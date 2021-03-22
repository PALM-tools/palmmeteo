import unittest


class TestCore(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_boolean_condition(self):
        self.assertTrue(True)

    def test_equality(self):
        self.assertEqual('one', 'one')
