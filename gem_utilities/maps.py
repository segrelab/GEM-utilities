import os
import time

import requests


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
