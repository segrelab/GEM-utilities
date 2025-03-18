import unittest

import cobra

from gem_utilities.media import clean_media


class TestCleanMedia(unittest.TestCase):
    def test_clean_media(self):
        # Load the textbook model
        model = cobra.io.load_model("textbook")

        # Define a media dictionary with a fake exchange reaction
        media = {
            "EX_glc__D_e": -10.0,  # valid exchange reaction
            "EX_o2_e": -20.0,  # valid exchange reaction
            "EX_foo_e": -5.0,  # fake exchange reaction
        }

        # Expected output should not include the fake exchange reaction
        expected_clean_media = {"EX_glc__D_e": -10.0, "EX_o2_e": -20.0}

        # Run the clean_media function
        clean_medium = clean_media(model, media)

        # Check if the cleaned media matches the expected output
        self.assertEqual(clean_medium, expected_clean_media)


if __name__ == "__main__":
    unittest.main()
