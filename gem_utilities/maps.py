import os
from typing import Dict, Iterable, List, Tuple

import cobra
import matplotlib.pyplot as plt
import seaborn as sns


class Mapper:
    """
    Make KEGG pathway maps incorporating data sourced from COBRApy models.

    Attributes
    ==========
    model : cobra.Model
        The COBRApy model containing the metabolic network.

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
        overwrite_output: bool = True,
        name_files: bool = False,
        categorize_files: bool = False,
    ) -> None:
        """
        Parameters
        ==========
        model : cobra.Model
            The COBRApy model containing the metabolic network.

        overwrite_output : bool, True
            If True, methods in this class overwrite existing output files.

        name_files : bool, False
            Include the pathway name along with the number in output map file names.

        categorize_files : bool, False
            Categorize output files by pathway map within subdirectories corresponding to the BRITE
            hierarchy of maps (see https://www.genome.jp/brite/br08901).
        """
        self.model = model
        self.overwrite_output = overwrite_output
        self.name_files = name_files
        self.categorize_files = categorize_files
        self.pathway_categorization = (
            self._categorize_pathways() if categorize_files else None
        )

    def map_reactions(
        self,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        color_hexcode: str = "#2ca02c",
        draw_maps_lacking_reactions: bool = False,
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting reactions present in the COBRA model.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            List of pathway numbers to draw. The default of None draws all available pathway maps.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions present in the model.

        draw_maps_lacking_reactions : bool, False
            If False, by default, only draw maps containing any of the reactions in the model.
            If True, draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select reactions and 'draw_maps_lacking_reactions' was
            False.
        """
        # Retrieve the IDs of all reactions in the model.
        reaction_ids = [reaction.id for reaction in self.model.reactions]

        drawn = self._map_reactions_fixed_colors(
            reaction_ids,
            output_dir,
            pathway_numbers=pathway_numbers,
            color_hexcode=color_hexcode,
            draw_maps_lacking_reactions=draw_maps_lacking_reactions,
        )
        count = sum(drawn.values()) if drawn else 0
        print(f"Number of maps drawn: {count}")

        return drawn

    def _map_reactions_fixed_colors(
        self,
        reaction_ids: Iterable[str],
        output_dir: str,
        pathway_numbers: List[str] = None,
        color_hexcode: str = "#2ca02c",
        draw_maps_lacking_reactions: bool = False,
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting reactions in a single color provided by a hex code.

        Parameters
        ==========
        reaction_ids : Iterable[str]
            Reaction IDs to be highlighted in the maps.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            List of pathway numbers to draw. The default of None draws all available pathway maps.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions present in the model.

        draw_maps_lacking_reactions : bool, False
            If False, by default, only draw maps containing any of the reactions in the model.
            If True, draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select reactions and 'draw_maps_lacking_reactions' was
            False.
        """
        # Find the numeric IDs of the maps to draw.
        pathway_numbers = self._find_maps(
            output_dir, "reactions", patterns=pathway_numbers
        )

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Draw maps.
        drawn: Dict[str, bool] = {}
        for pathway_number in pathway_numbers:
            drawn[pathway_number] = self._draw_map_reactions_single_color(
                pathway_number,
                reaction_ids,
                color_hexcode,
                output_dir,
                draw_maps_lacking_reactions,
            )

        return drawn

    def _draw_map_reactions_single_color(
        self,
        pathway_number: str,
        reaction_ids: Iterable[str],
        color_hexcode: str,
        output_dir: str,
        draw_map_lacking_reactions: bool = False,
    ) -> bool:
        """
        Draw a single pathway map, highlighting reactions in a single color.

        Parameters
        ==========
        pathway_number : str
            Pathway number to draw.

        reaction_ids : Iterable[str]
            Reaction IDs to be highlighted in the map.

        color_hexcode : str
            Color to use for highlighting reactions.

        output_dir : str
            Path to the output directory in which the pathway map PDF file is drawn.

        draw_map_lacking_reactions : bool
            If False, only draw maps containing any of the reactions in the model.
            If True, draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        bool
            True if the map was drawn, False otherwise.
        """
        # Placeholder for actual drawing logic
        # Here you would integrate with a library to draw the KEGG pathway maps
        # For now, we just simulate the drawing process
        print(f"Drawing map for pathway {pathway_number} with color {color_hexcode}")
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
            List of pathway numbers to draw. The default of None draws all available pathway maps.

        Returns
        =======
        List[str]
            List of pathway numbers to draw.
        """
        # Placeholder for actual logic to find maps
        # For now, we just return a list of dummy pathway numbers
        return patterns if patterns else ["00010", "00020", "00030"]

    def _categorize_pathways(self) -> dict:
        """
        Categorize pathways based on the BRITE hierarchy.

        Returns
        =======
        dict
            Dictionary mapping pathway numbers to categories.
        """
        # Placeholder for actual categorization logic
        # For now, we just return an empty dictionary
        return {}
