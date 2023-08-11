import unittest
from gem_utilities import helloworld
class MyTestCase(unittest.TestCase):
    def test_default_greeting_set(self):
        greeter = helloworld.Greeter()
        self.assertEqual(greeter.message, 'Hello world!')
if __name__ == '__main__':
    unittest.main()