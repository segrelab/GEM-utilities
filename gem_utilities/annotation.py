import time

import cobra
import requests


# Function to get KO numbers for a KEGG reaction ID using the KEGG API
def get_ko_for_kegg_reaction(kegg_reaction_id):
    # Remove 'R' prefix if present in a format like "R00001"
    if kegg_reaction_id.startswith("R"):
        kegg_id = kegg_reaction_id
    else:
        kegg_id = kegg_reaction_id

    # KEGG API URL
    url = f"http://rest.kegg.jp/link/ko/{kegg_id}"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            ko_list = []
            for line in response.text.strip().split("\n"):
                if line:  # Skip empty lines
                    parts = line.split("\t")
                    if len(parts) > 1:
                        ko = parts[1]
                        ko_list.append(ko)
            return ko_list
        else:
            return []
    except Exception as e:
        print(f"Error fetching KO for {kegg_id}: {e}")
        return []
    finally:
        # Be respectful of the KEGG API by adding a small delay
        time.sleep(0.5)


def add_kos_to_model(model: cobra.Model, verbose: bool = False):
    # Make a copy of the model to avoid modifying the original model
    working_model = model.copy()
    # List to store reactions without KEGG annotations
    reactions_without_kegg = []
    # Dictionary to cache KEGG to KO mappings
    # Cache to avoid repeated API calls
    kegg_to_ko_mapping = {}

    for reaction in working_model.reactions:
        kegg_ids = []

        # Check different possible annotation keys for KEGG reactions
        for key in ["kegg.reaction", "kegg", "kegg_reaction"]:
            if key in reaction.annotation:
                # Handle both string and list annotations
                if isinstance(reaction.annotation[key], str):
                    kegg_ids.append(reaction.annotation[key])
                else:
                    kegg_ids.extend(reaction.annotation[key])

        # Get KO numbers for the KEGG reactions you found
        if kegg_ids:
            ko_numbers = []
            for kegg_id in kegg_ids:
                # Use cached result if available
                if kegg_id in kegg_to_ko_mapping:
                    ko_numbers.extend(kegg_to_ko_mapping[kegg_id])
                else:
                    # Get KO numbers from KEGG API
                    kos = get_ko_for_kegg_reaction(kegg_id)
                    kegg_to_ko_mapping[kegg_id] = kos
                    ko_numbers.extend(kos)

            # Add unique KO numbers to reaction annotation
            if ko_numbers:
                reaction.annotation["ko"] = list(set(ko_numbers))
                print(
                    f"Added KO numbers for {reaction.id}: {reaction.annotation['ko']}"
                )
        else:
            reactions_without_kegg.append(reaction.id)

    if verbose:
        print(f"Reactions without KEGG annotations: {reactions_without_kegg}")

    return working_model
