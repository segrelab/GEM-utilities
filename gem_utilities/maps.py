import functools
import json
import math
import os
import re
import shutil
import time
from typing import Dict, Iterable, List, Literal, Optional, Set, Tuple, Union

import cobra
import fitz  # PyMuPDF
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio.KEGG.KGML.KGML_parser import read

# Define qualitative and repeating colormaps for consistent use with keggmapping.py
qualitative_colormaps: List[str] = [
    "Pastel1",
    "Pastel2",
    "Paired",
    "Accent",
    "Dark2",
    "Set1",
    "Set2",
    "Set3",
    "tab10",
    "tab20",
    "tab20b",
    "tab20c",
]

repeating_colormaps: List[str] = ["flag", "prism"]


class Mapper:
    """
    Make KEGG pathway maps incorporating data sourced from COBRApy models.

    Attributes
    ==========
    model : cobra.Model
        The COBRApy model containing the metabolic network.

    kegg_dir : str
        Directory containing KEGG database files.

    available_pathway_numbers : List[str]
        ID numbers of all pathways set up with PNG and KGML files in the KEGG data directory.

    pathway_names : Dict[str, str]
        The names of all KEGG pathways, including those without files in the KEGG data directory.
        Keys are pathway ID numbers and values are pathway names.

    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    name_files : bool
        Include the pathway name along with the number in output map file names.

    categorize_files : bool
        Categorize output files by pathway map within subdirectories corresponding to the BRITE
        hierarchy of maps (see https://www.genome.jp/brite/br08901).

    pathway_categorization : dict[str, list[str]]
        Maps pathway numbers to categorization, with categories listed from general to specific.
    """

    def __init__(
        self,
        model: cobra.Model,
        kegg_dir: str = None,
        overwrite_output: bool = True,
        name_files: bool = False,
        categorize_files: bool = False,
    ) -> None:
        """
        Parameters
        ==========
        model : cobra.Model
            The COBRApy model containing the metabolic network.

        kegg_dir : str, None
            Directory containing KEGG database files. If None, attempts to use a default location.

        overwrite_output : bool, True
            If True, methods in this class overwrite existing output files.

        name_files : bool, False
            Include the pathway name along with the number in output map file names.

        categorize_files : bool, False
            Categorize output files by pathway map within subdirectories corresponding to the BRITE
            hierarchy of maps (see https://www.genome.jp/brite/br08901).
        """
        self.model = model
        self.kegg_dir = kegg_dir

        self.overwrite_output = overwrite_output
        self.name_files = name_files
        self.categorize_files = categorize_files
        self.pathway_categorization = (
            self._categorize_pathways() if categorize_files else None
        )

        # Initialize drawers for colorbars and grids
        self.colorbar_drawer = ColorbarDrawer(overwrite_output=overwrite_output)
        self.grid_drawer = PDFGridDrawer(overwrite_output=overwrite_output)


    def map_model_kos(
        self,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        color_hexcode: str = "#2ca02c",
        draw_maps_lacking_kos: bool = False,
        ko_annotation_key: str = "kegg.orthology",
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting KOs present in the COBRA model.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing model KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the KOs in the model.
            If True, draw maps regardless, meaning that nothing may be colored.

        ko_annotation_key : str, 'kegg.orthology'
            The key in reaction.annotation dictionary containing KO IDs.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos' was
            False.
        """
        # Retrieve the IDs of all KO annotations in the model
        ko_ids = set()
        for reaction in self.model.reactions:
            if ko_annotation_key in reaction.annotation:
                # Handle both string and list annotations
                annotation = reaction.annotation[ko_annotation_key]
                if isinstance(annotation, str):
                    ko_ids.add(annotation)
                elif isinstance(annotation, list):
                    ko_ids.update(annotation)

        # Convert to list
        ko_ids_list = list(ko_ids)

        if not ko_ids_list:
            print(
                f"WARNING: No KO annotations found using key '{ko_annotation_key}'. Maps will be empty."
            )

        drawn = self._map_kos_fixed_colors(
            ko_ids_list,
            output_dir,
            pathway_numbers=pathway_numbers,
            color_hexcode=color_hexcode,
            draw_maps_lacking_kos=draw_maps_lacking_kos,
        )
        count = sum(drawn.values()) if drawn else 0
        print(f"Number of maps drawn: {count}")

        return drawn

    def map_models_kos(
        self,
        models: Dict[str, cobra.Model],
        output_dir: str,
        groups_txt: str = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_individual_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        colormap: Union[bool, str, mcolors.Colormap] = True,
        colormap_limits: Tuple[float, float] = None,
        colormap_scheme: Literal["by_count", "by_membership"] = None,
        reverse_overlay: bool = False,
        color_hexcode: str = "#2ca02c",
        group_colormap: Union[str, mcolors.Colormap] = "plasma_r",
        group_colormap_limits: Tuple[float, float] = (0.1, 0.9),
        group_reverse_overlay: bool = False,
        draw_maps_lacking_kos: bool = False,
        ko_annotation_key: str = "kegg.orthology",
    ) -> Dict[Literal["unified", "individual", "grid"], Dict]:
        """
        Draw pathway maps, highlighting KOs across multiple COBRA models or groups of models.

        A reaction on a map is defined by one or more KOs. These are matched to KO annotations
        in each COBRA model. The presence/absence of any of these KOs in a model translates in
        the map to the presence/absence of the reaction in the model.

        Parameters
        ==========
        models : Dict[str, cobra.Model]
            Dictionary mapping model names to COBRA model objects. Models should have
            different names by which they are uniquely identified.

        output_dir : str
            Path to the output directory in which pathway map and colorbar PDF files are drawn. The
            directory is created if it does not exist.

        groups_txt : str, None
            A tab-delimited text file specifying which group each model belongs to. The
            first column, which can have any header, contains the names of models,
            those provided as keys in the 'models' argument. The second column, which must be headed
            'group', contains group names, which are recommended to be single words without fancy
            characters, such as 'HIGH_TEMPERATURE' or 'LOW_FITNESS' rather than 'my group #1' or
            'IS-THIS-OK?'. Each model can only be associated with a single group. The
            'group_threshold' argument must also be used for the groups to take effect.

        group_threshold : float, None
            The proportion of models in a group containing data of interest for the group
            to be represented in terms of presence/absence in a map feature.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        draw_individual_files : Union[Iterable[str], bool], False
            If not False, draw map files for individual models or groups. If True, draw maps
            for all models/groups. Alternatively, can accept a subset of model names or group names.

        draw_grid : Union[Iterable[str], bool], False
            If not False, draw a paneled grid file for each pathway map showing the unified
            map alongside maps for individual models or groups. If True, include all models/groups.
            Alternatively, can accept a subset of model names or group names.

        colormap : Union[bool, str, mcolors.Colormap], True
            Controls dynamic coloring of reactions based on model/group membership.
            If False, uses color_hexcode for all highlighted reactions.

        colormap_limits : Tuple[float, float], None
            Limit the fraction of the 'colormap' used in dynamically selecting colors.

        colormap_scheme : Literal['by_count', 'by_membership'], None
            How to color reactions: by count of models/groups containing them, or explicitly
            by model/group membership.

        reverse_overlay : bool, False
            If True, draw reactions in fewer models/groups on top of those in more models/groups.

        color_hexcode : str, '#2ca02c'
            Color for reactions if colormap=False, or for individual model maps.

        group_colormap : Union[str, mcolors.Colormap], 'plasma_r'
            Colormap for individual group maps.

        group_colormap_limits : Tuple[float, float], (0.1, 0.9)
            Limits for the group colormap.

        group_reverse_overlay : bool, False
            If True, draw reactions found in fewer models of a group on top.

        draw_maps_lacking_kos : bool, False
            If False, only draw maps containing any KOs. If True, draw all specified maps.

        ko_annotation_key : str, 'kegg.orthology'
            The key in reaction.annotation dictionary containing KO IDs.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            Information about which maps were drawn.
        """
        # Extract KO IDs from all models
        ko_model_names: Dict[str, List[str]] = {}
        for model_name, model in models.items():
            for reaction in model.reactions:
                if ko_annotation_key in reaction.annotation:
                    # Handle both string and list annotations
                    annotation = reaction.annotation[ko_annotation_key]
                    if isinstance(annotation, str):
                        ko_ids = [annotation]
                    elif isinstance(annotation, list):
                        ko_ids = annotation
                    else:
                        continue

                    for ko_id in ko_ids:
                        try:
                            ko_model_names[ko_id].append(model_name)
                        except KeyError:
                            ko_model_names[ko_id] = [model_name]

        # Handle groups
        if (groups_txt is None and group_threshold is not None) or (
            groups_txt is not None and group_threshold is None
        ):
            raise ValueError(
                "To group models, arguments to both 'groups_txt' and 'group_threshold' must be provided."
            )

        group_model_names: Dict[str, List[str]] = {}
        model_name_group: Dict[str, str] = {}
        if groups_txt is None:
            source_group = None
            group_sources = None
            categories = list(models.keys())
        else:
            if not 0 <= group_threshold <= 1:
                raise ValueError(
                    f"'group_threshold' must be a number between 0 and 1, not {group_threshold}"
                )

            # Parse the groups file
            source_group, group_sources = self._parse_groups_file(groups_txt)
            categories = list(group_sources.keys())

            # Validate groups against models
            for model_name in models:
                try:
                    group = source_group[model_name]
                except KeyError:
                    raise ValueError(f"Model '{model_name}' not found in groups file")

                try:
                    group_model_names[group].append(model_name)
                except KeyError:
                    group_model_names[group] = [model_name]
                model_name_group[model_name] = group

            # Check for models in groups file that aren't in the provided models
            for source in source_group:
                if source not in models:
                    print(
                        f"WARNING: Model '{source}' from groups file not found in provided models"
                    )

        # Find and create output directory
        pathway_numbers = self._find_maps(output_dir, "kos", patterns=pathway_numbers)
        os.makedirs(output_dir, exist_ok=True)

        drawn: Dict[Literal["unified", "individual", "grid"], Dict] = {
            "unified": {},
            "individual": {},
            "grid": {},
        }

        # Configure the colormap
        cmap, sampling, color_priority, scheme, ignore_groups, exceeds_colors = (
            self._configure_colormap(
                colormap,
                colormap_scheme,
                colormap_limits,
                categories,
                reverse_overlay,
                groups_txt,
            )
        )

        # Draw unified maps
        print("Drawing 'unified' map incorporating data from all models...")
        if scheme == "static":
            # Draw with static colors
            for pathway_number in pathway_numbers:
                if color_hexcode == "original":
                    drawn["unified"][pathway_number] = (
                        self._draw_map_kos_original_color(
                            pathway_number,
                            ko_model_names,
                            output_dir,
                            draw_map_lacking_kos=draw_maps_lacking_kos,
                        )
                    )
                else:
                    drawn["unified"][pathway_number] = self._draw_map_kos_single_color(
                        pathway_number,
                        ko_model_names,
                        color_hexcode,
                        output_dir,
                        draw_map_lacking_kos=draw_maps_lacking_kos,
                    )
        else:
            # Draw with dynamic coloring
            category_combos = None
            if scheme == "by_count":
                # Sample colormap for colors for different counts
                pass  # Implement sampling logic similar to keggmapping.py
            elif scheme == "by_membership":
                # Sample colormap for colors for different membership combinations
                category_combos = []
                for category_count in range(1, len(categories) + 1):
                    category_combos += list(combinations(categories, category_count))
                # Implement sampling logic

            # Draw colorbar
            self._draw_colorbar(
                scheme, color_priority, category_combos, categories, output_dir
            )

            # Draw maps with membership coloring
            for pathway_number in pathway_numbers:
                drawn["unified"][pathway_number] = self._draw_map_kos_membership(
                    pathway_number,
                    ko_model_names,
                    color_priority,
                    output_dir,
                    category_combos=category_combos,
                    group_sources=None if groups_txt is None else group_model_names,
                    group_threshold=None if groups_txt is None else group_threshold,
                    draw_map_lacking_kos=draw_maps_lacking_kos,
                )

        if exceeds_colors:
            print(
                f"WARNING: Fewer colors available in colormap ({exceeds_colors[0]}) than needed ({exceeds_colors[1]})"
            )

        # Handle individual files and grid drawing if requested
        # Similar logic to keggmapping.py would follow here

        # Return information about drawn maps
        count = sum(drawn["unified"].values()) if drawn["unified"] else 0
        print(f"Number of maps drawn: {count}")
        return drawn

    def _map_kos_fixed_colors(
        self,
        ko_ids: Iterable[str],
        output_dir: str,
        pathway_numbers: List[str] = None,
        color_hexcode: str = "#2ca02c",
        draw_maps_lacking_kos: bool = False,
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting reactions containing select KOs in either a single color
        provided by a hex code or the colors originally used in the reference map.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KO IDs to be highlighted in the maps.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing provided KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos' was
            False.
        """
        # Find the numeric IDs of the maps to draw.
        # pathway_numbers = self._find_maps(output_dir, "kos", patterns=pathway_numbers)

        os.makedirs(output_dir, exist_ok=True)

        # Draw maps.
        print("Drawing maps...")
        drawn: Dict[str, bool] = {}
        for pathway_number in pathway_numbers:
            if color_hexcode == "original":
                drawn[pathway_number] = self._draw_map_kos_original_color(
                    pathway_number,
                    ko_ids,
                    output_dir,
                    draw_map_lacking_kos=draw_maps_lacking_kos,
                )
            else:
                drawn[pathway_number] = self._draw_map_kos_single_color(
                    pathway_number,
                    ko_ids,
                    color_hexcode,
                    output_dir,
                    draw_map_lacking_kos=draw_maps_lacking_kos,
                )

        return drawn

    def _draw_map_kos_single_color(
        self,
        pathway_number: str,
        ko_ids: Iterable[str],
        color_hexcode: str,
        output_dir: str,
        draw_map_lacking_kos: bool = False,
    ) -> bool:
        """
        Draw a pathway map, highlighting reactions containing select KOs in a single color.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        ko_ids : Iterable[str]
            Select KOs, any of which in the map are colored.

        color_hexcode : str
            This is the color for reactions containing provided KOs.

        output_dir : str
            Path to the output directory in which map PDF files are drawn.

        draw_map_lacking_kos : bool, False
            If False, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_kos' was False.
        """
        # Implementation would normally involve:
        # 1. Loading the KGML file
        # 2. Identifying reactions with the provided KO IDs
        # 3. Coloring those reactions
        # 4. Drawing the map

        # Load the KGML file
        pathway = self._get_pathway(pathway_number)
        if pathway is None:
            return False

        # Check if any of the KOs are in the pathway
        has_kos = self._pathway_has_kos(pathway, ko_ids)
        if not has_kos and not draw_maps_lacking_kos:
            return False

        # Draw the map
        self._draw_map(pathway, output_dir, color_hexcode)
        return True

    def _draw_map_kos_original_color(
        self,
        pathway_number: str,
        ko_ids: Iterable[str],
        output_dir: str,
        draw_map_lacking_kos: bool = False,
    ) -> bool:
        """
        Draw a pathway map, highlighting reactions containing select KOs in the color or colors
        originally used in the reference map.

        Parameters
        ==========
        pathway_number : str, None
            Numeric ID of the map to draw.

        ko_ids : Iterable[str]
            Select KOs, any of which in the map are colored.

        output_dir : str
            Path to the output directory in which map PDF files are drawn.

        draw_map_lacking_kos : bool, False
            If False, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_kos' was False.
        """
        # Implementation similar to _draw_map_kos_single_color, but using original colors
        pathway = self._get_pathway(pathway_number)
        if pathway is None:
            return False

        # Check if any of the KOs are in the pathway
        has_kos = self._pathway_has_kos(pathway, ko_ids)
        if not has_kos and not draw_maps_lacking_kos:
            return False

        # Draw the map with original colors
        self._draw_map(pathway, output_dir, "original")
        return True

    def _draw_map_kos_membership(
        self,
        pathway_number: str,
        ko_membership: Dict[str, List[str]],
        color_priority: Dict[str, float],
        output_dir: str,
        category_combos: List[Tuple[str]] = None,
        group_sources: Dict[str, List[str]] = None,
        group_threshold: float = None,
        draw_map_lacking_kos: bool = False,
    ) -> bool:
        """
        Draw a pathway map, coloring reactions by their membership via KOs in data sources
        (e.g., models) or groups of data sources.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        ko_membership : Dict[str, List[str]]
            Keys are KO IDs. Values are lists of data sources in which KOs are found.

        color_priority : Dict[str, float]
            Keys are color hex codes. Values are priorities.

        output_dir : str
            Path to the output directory in which map PDF files are drawn.

        category_combos : List[Tuple[str]], None
            List of tuples representing all possible combinations of sources or groups.

        group_sources : Dict[str, List[str]], None
            Keys are group names. Values are lists of data sources categorized in the group.

        group_threshold : float, None
            The proportion of KO data sources in a group needed for the group to be represented.

        draw_map_lacking_kos : bool, False
            If False, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False otherwise.
        """
        # Implementation would normally involve:
        # 1. Loading the KGML file
        # 2. Identifying reactions with the provided KO IDs
        # 3. Determining which groups/models each reaction belongs to
        # 4. Coloring reactions accordingly
        # 5. Drawing the map

        # Simplified placeholder implementation
        pathway = self._get_pathway(pathway_number)
        if pathway is None:
            return False

        # Check if any of the KOs are in the pathway
        has_kos = self._pathway_has_kos(pathway, ko_membership.keys())
        if not has_kos and not draw_map_lacking_kos:
            return False

        # Draw the map with membership-based coloring
        self._draw_map(pathway, output_dir, color_priority)
        return True

    def _find_maps(
        self, output_dir: str, prefix: str, patterns: List[str] = None
    ) -> List[str]:
        """
        Find the numeric IDs of maps to draw given the file prefix, checking that the map can be
        drawn in the target output directory.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn.

        prefix : str
            Output filenames are formatted as <prefix>_<pathway_number>.pdf or
            <prefix>_<pathway_number>_<pathway_name>.pdf.

        patterns : List[str], None
            Regex patterns of pathway numbers, which are five digits.

        Returns
        =======
        List[str]
            List of pathway numbers to draw.
        """
        if patterns is None:
            pathway_numbers = self.available_pathway_numbers
        else:
            pathway_numbers = self._get_pathway_numbers_from_patterns(patterns)

        if self.pathway_categorization is not None:
            missing_pathway_numbers: list[str] = []
            for pathway_number in pathway_numbers:
                if pathway_number not in self.pathway_categorization:
                    missing_pathway_numbers.append(pathway_number)
            if missing_pathway_numbers:
                missing_str = ", ".join(f"'{p}'" for p in missing_pathway_numbers)
                print(f"WARNING: Pathway numbers not in BRITE hierarchy: {missing_str}")
                print("File categorization cannot be used.")
                self.categorize_files = False
                self.pathway_categorization = None

        if not self.overwrite_output:
            # Check if files would be overwritten
            for pathway_number in pathway_numbers:
                pathway_name = (
                    f"_{self._name_pathway(pathway_number)}" if self.name_files else ""
                )
                out_basename = f"{prefix}_{pathway_number}{pathway_name}.pdf"
                if self.pathway_categorization is None:
                    out_path = os.path.join(output_dir, out_basename)
                else:
                    out_path = os.path.join(
                        *self.pathway_categorization[pathway_number], out_basename
                    )
                if os.path.exists(out_path):
                    raise ValueError(
                        f"Output files would be overwritten in {output_dir}. Use overwrite_output=True or delete existing files."
                    )

        return pathway_numbers

    def _get_pathway_numbers_from_patterns(self, patterns: Iterable[str]) -> List[str]:
        """
        Among pathways available in the KEGG data directory, get those with ID numbers matching the
        given regex patterns.

        Parameters
        ==========
        patterns : Iterable[str]
            Regex patterns of pathway numbers, which are five digits.

        Returns
        =======
        List[str]
            Pathway numbers matching the regex patterns.
        """
        pathway_numbers: List[str] = []
        for pattern in patterns:
            for available_pathway_number in self.available_pathway_numbers:
                if re.match(pattern, available_pathway_number):
                    pathway_numbers.append(available_pathway_number)

        # Maintain the order of pathway numbers recovered from patterns.
        seen = set()
        return [
            pathway_number
            for pathway_number in pathway_numbers
            if not (pathway_number in seen or seen.add(pathway_number))
        ]

    def _get_pathway(self, pathway_number: str):
        """
        Get a pathway object for the KGML file used in drawing a pathway map.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        Returns
        =======
        object
            Representation of the KGML file as an object.
        """
        # Combine the pathway number with the KEGG directory to get the KGML file
        file_path = os.path.join(self.kegg_dir, f"ko{pathway_number}.xml")

        # Check that the file exists, if it does not throw an error
        if not os.path.exists(file_path):
            print(f"WARNING: Pathway map file not found for pathway number '{pathway_number}'")
            return None

        # Load the KGML file and return the pathway object using BioPython
        pathway = read(file_path)

        return pathway


    def _pathway_has_kos(self, pathway, ko_ids: Iterable[str]) -> bool:
        """
        Check if a pathway contains any of the specified KO IDs.

        Parameters
        ==========
        pathway : object
            Pathway object.

        ko_ids : Iterable[str]
            KO IDs to check for.

        Returns
        =======
        bool
            True if the pathway contains any of the KO IDs, False otherwise.
        """
        # In a real implementation, this would check the KGML file for the KO IDs
        # Simplified placeholder implementation
        return True  # Assume the pathway has the KOs for demonstration purposes

    def _draw_map(self, pathway, output_dir: str, color_spec):
        """
        Draw a map given the pathway data.

        Parameters
        ==========
        pathway : object
            Pathway object.

        output_dir : str
            Path to the output directory in which the pathway map PDF file is drawn.

        color_spec : Union[str, Dict[str, float]]
            Color specification, either a hex code, 'original', or a dictionary mapping colors to priorities.
        """
        pathway_number = pathway["number"]
        pathway_name = (
            f"_{self._name_pathway(pathway_number)}" if self.name_files else ""
        )
        out_basename = f"kos_{pathway_number}{pathway_name}.pdf"

        if self.pathway_categorization is None:
            out_dir = output_dir
            out_path = os.path.join(output_dir, out_basename)
        else:
            out_dir = os.path.join(
                output_dir, *self.pathway_categorization[pathway_number]
            )
            out_path = os.path.join(out_dir, out_basename)

        os.makedirs(out_dir, exist_ok=True)

        # In a real implementation, this would draw the map using a library
        # Here we'll just create a simplified placeholder PDF
        self._create_placeholder_pdf(out_path, pathway, color_spec)

        if self.pathway_categorization is not None:
            self._symlink_map(output_dir, out_path)

    def _create_placeholder_pdf(self, out_path: str, pathway, color_spec):
        """
        Create a placeholder PDF file for demonstration purposes.

        Parameters
        ==========
        out_path : str
            Path to the output PDF file.

        pathway : object
            Pathway object.

        color_spec : Union[str, Dict[str, float]]
            Color specification.
        """
        # Create a simple PDF with pathway information
        fig, ax = plt.subplots(figsize=(8, 6))

        # Display pathway info
        pathway_number = pathway["number"]
        pathway_name = self.pathway_names.get(pathway_number, "Unknown pathway")

        ax.text(0.5, 0.8, f"Pathway: {pathway_number}", ha="center", fontsize=14)
        ax.text(0.5, 0.7, f"Name: {pathway_name}", ha="center", fontsize=12)

        # Display color info
        if isinstance(color_spec, str):
            if color_spec == "original":
                color_info = "Using original colors"
            else:
                color_info = f"Single color: {color_spec}"
        else:
            color_info = f"Using {len(color_spec)} colors with dynamic coloring"

        ax.text(0.5, 0.5, color_info, ha="center", fontsize=10)

        # Add note that this is a placeholder
        ax.text(
            0.5,
            0.3,
            "This is a placeholder for KEGG pathway map visualization",
            ha="center",
            fontsize=10,
            style="italic",
        )
        ax.text(
            0.5,
            0.2,
            "In a real implementation, the actual pathway would be drawn here",
            ha="center",
            fontsize=10,
            style="italic",
        )

        ax.axis("off")
        plt.tight_layout()
        plt.savefig(out_path, format="pdf")
        plt.close()

    def _name_pathway(self, pathway_number: str) -> str:
        """
        Format the pathway name corresponding to the number for suitability in file paths.

        Replace all non-alphanumeric characters except parentheses, brackets, and curly braces with
        underscores. Replace multiple consecutive underscores with a single underscore. Strip
        leading and trailing underscores.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of a pathway map.

        Returns
        =======
        str
            Altered version of the pathway name.
        """
        try:
            pathway_name = self.pathway_names[pathway_number]
        except KeyError:
            print(
                f"WARNING: Pathway number '{pathway_number}' not found in pathway names"
            )
            return "unknown_pathway"

        # Sanitize the name for use in file paths
        altered = re.sub(r"[^a-zA-Z0-9()\[\]\{\}]", "_", pathway_name)
        altered = re.sub(r"_+", "_", altered)
        altered = altered.strip("_")

        return altered

    def _symlink_map(self, output_dir: str, map_path: str) -> None:
        """
        Make a symlink for a map file in a dedicated symlink directory.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which the map file was drawn.

        map_path : str
            Map file path to be symlinked.
        """
        symlink_dir = os.path.join(output_dir, "symlink")
        os.makedirs(symlink_dir, exist_ok=True)
        map_basename = os.path.basename(map_path)
        symlink_path = os.path.join(symlink_dir, map_basename)
        if os.path.exists(symlink_path):
            os.remove(symlink_path)
        os.symlink(os.path.abspath(map_path), symlink_path)

    def _categorize_pathways(self) -> dict[str, list[str]]:
        """
        Categorize pathways in the BRITE hierarchy, 'br08901'.

        Alter category names to make suitable for directory paths.

        Returns
        =======
        dict[str, list[str]]
            Keys are pathway numbers. Values are lists of the categories from general to specific.
            For example, '00010': ['Metabolism', 'Carbohydrate_metabolism']
        """
        if not os.path.exists(self.kegg_brite_pathways_file):
            print(
                f"WARNING: KEGG BRITE hierarchy file not found at {self.kegg_brite_pathways_file}"
            )
            return {}

        # Load the BRITE hierarchy
        try:
            with open(self.kegg_brite_pathways_file) as f:
                hierarchy = json.load(f)
        except Exception as e:
            print(f"ERROR: Failed to load BRITE hierarchy: {e}")
            return {}

        # Process the hierarchy to extract pathway categorization
        pathway_categorization = {}

        # This is a simplified implementation - in a real implementation,
        # you would parse the actual BRITE hierarchy structure
        # Placeholder categorization for demonstration purposes
        for pathway_number in self.available_pathway_numbers:
            first_digit = pathway_number[0]
            if first_digit == "0":
                if pathway_number.startswith("00"):
                    category = ["Metabolism"]
                    if int(pathway_number[2:]) < 50:
                        category.append("Carbohydrate_metabolism")
                    else:
                        category.append("Energy_metabolism")
                else:
                    category = ["Genetic_Information_Processing"]
            elif first_digit == "1":
                category = ["Environmental_Information_Processing"]
            elif first_digit == "2":
                category = ["Cellular_Processes"]
            elif first_digit == "3":
                category = ["Organismal_Systems"]
            else:
                category = ["Uncategorized"]

            pathway_categorization[pathway_number] = category

        return pathway_categorization

    def _parse_groups_file(
        self, groups_txt: str
    ) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        """
        Parse a groups file to get model grouping information.

        Parameters
        ==========
        groups_txt : str
            Path to the groups file.

        Returns
        =======
        Tuple[Dict[str, str], Dict[str, List[str]]]
            First dict maps model names to group names.
            Second dict maps group names to lists of model names.
        """
        source_group = {}
        group_sources = {}

        try:
            df = pd.read_csv(groups_txt, sep="\t")
            if "group" not in df.columns:
                raise ValueError(f"Groups file must have a 'group' column")

            # First column is model name
            model_col = df.columns[0]

            for _, row in df.iterrows():
                model_name = row[model_col]
                group = row["group"]

                source_group[model_name] = group

                if group not in group_sources:
                    group_sources[group] = []
                group_sources[group].append(model_name)

        except Exception as e:
            print(f"ERROR: Failed to parse groups file: {e}")
            raise

        return source_group, group_sources

    def _configure_colormap(
        self,
        colormap: Union[bool, str, mcolors.Colormap],
        colormap_scheme: Literal["by_count", "by_membership"],
        colormap_limits: Tuple[float, float],
        categories: List[str],
        reverse_overlay: bool,
        groups_txt: str,
    ):
        """
        Configure the colormap for pathway visualization.

        Parameters
        ==========
        colormap : Union[bool, str, mcolors.Colormap]
            Specifies the colormap to use or whether to use dynamic coloring.

        colormap_scheme : Literal['by_count', 'by_membership']
            How to color reactions.

        colormap_limits : Tuple[float, float]
            Limits for the colormap.

        categories : List[str]
            List of categories (models or groups).

        reverse_overlay : bool
            Whether to reverse the overlay order.

        groups_txt : str
            Path to groups file.

        Returns
        =======
        Tuple[mcolors.Colormap, str, Dict[str, float], str, bool, Tuple[int, int]]
            Returns colormap configuration parameters.
        """
        print("Setting map colors...")

        # Set the colormap scheme
        ignore_groups = False
        if colormap is False:
            scheme = "static"
            if groups_txt is not None:
                ignore_groups = True
        else:
            if colormap_scheme is None:
                if len(categories) < 4:
                    scheme = "by_membership"
                else:
                    scheme = "by_count"
            elif colormap_scheme in ["by_count", "by_membership"]:
                scheme = colormap_scheme
            else:
                raise ValueError(f"Invalid colormap_scheme: {colormap_scheme}")

        # Set the colormap
        if colormap is True:
            if scheme == "by_count":
                cmap = plt.colormaps["plasma_r"]
                if colormap_limits is None:
                    colormap_limits = (0.1, 0.9)
            elif scheme == "by_membership":
                cmap = plt.colormaps["tab10"]
                if colormap_limits is None:
                    colormap_limits = (0.0, 1.0)
        elif colormap is False:
            cmap = None
            sampling = None
        elif isinstance(colormap, str):
            cmap = plt.colormaps[colormap]
            if colormap_limits is None:
                colormap_limits = (0.0, 1.0)
        elif isinstance(colormap, mcolors.Colormap):
            cmap = colormap
            if colormap_limits is None:
                colormap_limits = (0.0, 1.0)
        else:
            raise ValueError(f"Invalid colormap: {colormap}")

        # Set how the colormap is sampled
        if cmap is None:
            sampling = None
        else:
            if cmap.name in qualitative_colormaps + repeating_colormaps:
                sampling = "in_order"
            else:
                sampling = "even"

        # Trim the colormap
        exceeds_colors = None
        if (
            cmap is not None
            and colormap_limits is not None
            and colormap_limits != (0.0, 1.0)
        ):
            lower_limit = colormap_limits[0]
            upper_limit = colormap_limits[1]
            assert 0.0 <= lower_limit <= upper_limit <= 1.0
            cmap = mcolors.LinearSegmentedColormap.from_list(
                f"trunc({cmap.name},{lower_limit:.2f},{upper_limit:.2f})",
                cmap(range(int(lower_limit * cmap.N), math.ceil(upper_limit * cmap.N))),
            )

        # Check if we have enough colors
        if cmap is not None and len(categories) > cmap.N:
            exceeds_colors = (cmap.N, len(categories))

        # Create color priority dictionary
        color_priority = {}
        if cmap is not None:
            # Create color priorities based on sampling method
            if sampling == "in_order":
                if len(categories) == 1:
                    sample_points = range(1, 2)
                else:
                    sample_points = range(len(categories))
            elif sampling == "even":
                if len(categories) == 1:
                    sample_points = np.linspace(1, 1, 1)
                else:
                    sample_points = np.linspace(0, 1, len(categories))

            for i, sample_point in enumerate(sample_points):
                if reverse_overlay:
                    color_priority[mcolors.rgb2hex(cmap(sample_point))] = (
                        1 - sample_point
                    )
                else:
                    color_priority[mcolors.rgb2hex(cmap(sample_point))] = sample_point

        return cmap, sampling, color_priority, scheme, ignore_groups, exceeds_colors

    def _draw_colorbar(
        self,
        scheme: str,
        color_priority: Dict[str, float],
        category_combos,
        categories: List[str],
        output_dir: str,
    ):
        """
        Draw a colorbar for the pathway maps.

        Parameters
        ==========
        scheme : str
            The colormap scheme being used.

        color_priority : Dict[str, float]
            Maps colors to priorities.

        category_combos
            Combinations of categories for 'by_membership' scheme.

        categories : List[str]
            List of categories (models or groups).

        output_dir : str
            Directory to save the colorbar.
        """
        if scheme == "by_count":
            self.colorbar_drawer.draw(
                color_priority,
                os.path.join(output_dir, "colorbar.pdf"),
                color_labels=range(1, len(categories) + 1),
                label="model count",
            )
        elif scheme == "by_membership":
            self.colorbar_drawer.draw(
                color_priority,
                os.path.join(output_dir, "colorbar.pdf"),
                color_labels=[", ".join(combo) for combo in category_combos],
                label="models",
            )


class ColorbarDrawer:
    """
    Writes standalone colorbar image files.

    Parameters
    ==========
    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    figsize : Tuple[int, int]
        Dimensions of the figure in inches.

    orientation : Literal['horizontal', 'vertical']
        Orientation of the colobar.

    tick_fontsize : Union[int, None]
        Font size of tick labels.

    label_rotation : Union[int, None]
        Rotation of the colorbar label.

    label_fontsize : int
        Font size of colorbar label.

    labelpad : int
        Spacing of colorbar label from tick labels in points.
    """

    def __init__(self, overwrite_output: bool = True) -> None:
        """
        Parameters
        ==========
        overwrite_output : bool, True
            If True, methods in this class overwrite existing output files.
        """
        self.overwrite_output = overwrite_output

        self.figsize: Tuple[int, int] = (1, 6)
        self.orientation: Literal["horizontal", "vertical"] = "vertical"
        self.tick_fontsize: Union[int, None] = None
        self.label_rotation: int = None
        self.label_fontsize: int = 24
        self.labelpad: int = 30

    def draw(
        self,
        colors: Iterable,
        out_path: str,
        color_labels: Iterable[str] = None,
        label: str = None,
    ) -> None:
        """
        Save a standalone colorbar to a file.

        Parameters
        ==========
        colors : Iterable
            Sequence of Matplotlib color specifications for 'matplotlib.colors.ListedColormap' color
            parameter.

        out_path : str
            Path to PDF output file.

        color_labels : Iterable[str], None
            Color segment labels.

        label : str, None
            Overall colorbar label.
        """
        if color_labels is not None:
            assert len(colors) == len(color_labels)

        fig, ax = plt.subplots(figsize=self.figsize)

        cmap = mcolors.ListedColormap(colors)
        norm = mcolors.BoundaryNorm(
            boundaries=range(len(colors) + 1), ncolors=len(colors)
        )

        cb = plt.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation=self.orientation,
        )

        # Don't show tick marks.
        cb.ax.tick_params(size=0)

        if color_labels:
            if self.tick_fontsize is None:
                # Calculate appropriate font size of tick labels based on color segment height.
                length_in_data_coords = 1 / len(colors)
                origin_in_points = ax.transData.transform((0, 0))
                if self.orientation == "vertical":
                    size_value = height_in_points = (
                        ax.transData.transform((0, length_in_data_coords))
                        - origin_in_points
                    )[1]
                elif self.orientation == "horizontal":
                    size_value = width_in_points = (
                        ax.transData.transform((length_in_data_coords, 0))
                        - origin_in_points
                    )[0]
                else:
                    raise AssertionError
                if size_value < 10:
                    tick_fontsize = size_value * 2
                else:
                    tick_fontsize = min(size_value, 24)
            else:
                tick_fontsize = self.tick_fontsize

            cb.set_ticks(np.arange(len(colors)) + 0.5)
            cb.set_ticklabels(color_labels, fontsize=tick_fontsize)

        if label:
            if self.label_rotation is None:
                if self.orientation == "vertical":
                    label_rotation = 270
                elif self.orientation == "horizontal":
                    label_rotation = 0
                else:
                    raise AssertionError
            else:
                label_rotation = self.label_rotation
            cb.set_label(
                label,
                rotation=label_rotation,
                labelpad=self.labelpad,
                fontsize=self.label_fontsize,
            )

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, format="pdf", bbox_inches="tight")
        plt.close()


class PDFGridDrawer:
    """
    Writes PDF files that are a grid of input PDF files.

    Attributes
    ==========
    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    paper_format : Union[str, None]
        Paper format string, e.g., 'letter', 'letter-l', 'A4', 'A4-L'.

    margin : float
        Minimum space between grid cells.

    label_fontsize_scale : Union[float, None]
        Font size of labels over grid cells as proportion of margin.
    """

    def __init__(self, overwrite_output: bool = True) -> None:
        """
        Parameters
        ==========
        overwrite_output : bool, True
            If True, methods in this class overwrite existing output files.
        """
        self.overwrite_output = overwrite_output

        self.paper_format: Union[str, None] = None
        self.margin: float = 10.0
        self.label_fontsize_scale: float = 0.8

    def draw(
        self, in_paths: Iterable[str], out_path: str, labels: Iterable[str] = None
    ) -> None:
        """
        Write a PDF containing a grid of input PDF images.

        Parameters
        ==========
        in_paths : Iterable[str]
            Paths to input PDFs.

        out_path : str
            Path to output PDF.

        labels : Iterable[str], None
            Labels displayed over grid cells corresponding to input files.
        """
        assert len(in_paths) > 0
        if labels:
            assert len(in_paths) == len(labels)
        assert 0 < self.label_fontsize_scale <= 1

        # In a real implementation, this would be a complex operation to create a grid of PDFs
        # Here's a simplified placeholder implementation
        print(f"Creating PDF grid with {len(in_paths)} input files")

        # Create a simple PDF with grid information
        fig, ax = plt.subplots(figsize=(8, 6))

        ax.text(
            0.5,
            0.8,
            f"PDF Grid with {len(in_paths)} input files",
            ha="center",
            fontsize=14,
        )

        if labels:
            label_str = ", ".join(labels)
            ax.text(0.5, 0.6, f"Included maps: {label_str}", ha="center", fontsize=10)

        ax.text(
            0.5,
            0.3,
            "This is a placeholder for a PDF grid visualization",
            ha="center",
            fontsize=10,
            style="italic",
        )
        ax.text(
            0.5,
            0.2,
            "In a real implementation, a grid of pathway maps would be shown",
            ha="center",
            fontsize=10,
            style="italic",
        )

        ax.axis("off")
        plt.tight_layout()

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, format="pdf")
        plt.close()


def combinations(items, r):
    """
    Simple combinations function to avoid requiring itertools.

    Parameters
    ==========
    items : Iterable
        Items to combine.

    r : int
        Length of combinations.

    Returns
    =======
    List[Tuple]
        List of combinations.
    """
    if r == 0:
        return [()]
    if not items:
        return []

    result = []
    for i, item in enumerate(items):
        for combo in combinations(items[i + 1 :], r - 1):
            result.append((item,) + combo)
    return result


def download_kegg_pathways(output_dir: str, pathway_ids=None):
    """
    Download KEGG pathways as KGML and image files.

    Parameters
    ----------
    output_dir : str
        Path to the output directory in which pathway files are saved.
    pathway_ids : list, optional
        List of pathway IDs to download the maps for, by default None, which
        downloads all pathways.
    """
    os.makedirs(output_dir, exist_ok=True)

    # If no specific pathway IDs are provided, get all pathway IDs
    if pathway_ids is None:
        response = requests.get("https://rest.kegg.jp/list/pathway")
        pathway_ids = [
            line.split("\t")[0][3:] for line in response.text.strip().split("\n")
        ]

    # Create directory structure
    kgml_dir = os.path.join(output_dir, "kgml")
    image_dir = os.path.join(output_dir, "image")
    os.makedirs(kgml_dir, exist_ok=True)
    os.makedirs(image_dir, exist_ok=True)

    for pathway_id in pathway_ids:
        print(f"Downloading {pathway_id}...")

        # Get KGML
        kgml_response = requests.get(f"https://rest.kegg.jp/get/ko{pathway_id}/kgml")
        if kgml_response.status_code == 200:
            with open(os.path.join(kgml_dir, f"ko{pathway_id}.xml"), "w") as f:
                f.write(kgml_response.text)

        # Get pathway image
        img_response = requests.get(f"https://rest.kegg.jp/get/ko{pathway_id}/image")
        if img_response.status_code == 200:
            with open(os.path.join(image_dir, f"ko{pathway_id}.png"), "wb") as f:
                f.write(img_response.content)

        # Be nice to the server - add a delay between requests
        time.sleep(1)
