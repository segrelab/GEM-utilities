import os
import time
from typing import List

import requests
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import kegg_get


def map_ko_ids(
    pathway_id: str,
    ko_ids: List[str],
    color="#BFBFFF",
    kgml_folder=".",
    output_folder=".",
):
    """
    Render a KEGG pathway map with specified KO identifiers highlighted in the
    given color, save as a PDF and a PNG.

    Parameters
    ----------
    pathway_id : str
        Pathway ID (e.g. "ko00061")
    ko_ids : _type_, List[str]
        List of KEGG Ortholog IDs to highlight on the map
        (e.g. ["K00665", "K00668"])
    color : str, optional
        Hex code for color to highlight the KOs with, by default "#BFBFFF"
    kgml_folder : str, optional
        Path to look for and/or save KGML files, by default '.'
    output_folder : str, optional
        Path of folder to save resulting PDFs and PNGs, by default "."
    """
    # If the pathway's KGML file is not already downloaded, download it
    if not os.path.exists(os.path.join(kgml_folder, pathway_id + ".xml")):
        with open(os.path.join(kgml_folder, pathway_id + ".xml"), "w") as f:
            f.write(kegg_get(pathway_id, "kgml").read())

    # Get the path to the KGML file
    file_path = os.path.join(kgml_folder, pathway_id + ".xml")

    # Read the KGML file
    with open(file_path, "r") as f:
        pathway = KGML_parser.read(f)

    # Make all squares white, except for the ones with the KO identifiers
    for element in pathway.orthologs:
        for graphic in element.graphics:
            if graphic.name in ko_ids:
                graphic.bgcolor = color
            else:
                graphic.bgcolor = "#FFFFFF"

    # Render the map
    canvas = KGMLCanvas(pathway, import_imagemap=True)

    # Save the file as a PDF and convert it to a PNG
    canvas.draw(os.path.join(output_folder, pathway_id + ".pdf"))
    os.system("convert {0}.pdf {0}.png".format(os.path.join(output_folder, pathway_id)))


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
