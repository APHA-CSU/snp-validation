import unittest
from unittest.mock import Mock

import validator
from samples.sample import Sample
import entities.sequenced as sequenced
import entities.processed as processed
import entities.genome as genome
import utils.utils as utils

class ValidatorTests(unittest.TestCase):
    def test_validator(self):
        """ A simple introductory unit test that checks we can run through the
            validator without errors
        """
        # Mock I/O on validator
        validator.shutil.copytree = Mock()
        validator.os.makedirs = Mock()
        validator.os.path.exists = Mock(return_value=False)
        utils.utils.run = Mock()
        utils.utils.glob.glob = Mock(return_value=[''])

        validator.compare_snps.benchmark = Mock(return_value=[Mock(), Mock(), {'mock': Mock()}])

        # Mock Sample
        sequenced.from_results_dir = Mock(return_value=[Mock()])
        processed.from_list = Mock(return_value=[Mock()])
        
        mock_sample = Sample()

        genome.os.path.exists = Mock(return_value=True)
        mock_genome = genome.SimulatedGenome('mock-genome', '', '', '', '')

        mock_sample.simulate_genome = Mock(return_value=mock_genome)
        mock_sample.simulate_reads = Mock(return_value=['READ1', 'READ2'])

        # Run test
        validator.pipeline('./', './', [mock_sample])

    def test_names_consistent(self):
        """ Ensure names consistent works as expected """
        self.assertTrue(utils.names_consistent([], []))
        self.assertTrue(utils.names_consistent([mock_name("A")], [mock_name("A")]))
        self.assertTrue(utils.names_consistent(mock_names("AB"), mock_names("BA")))
        
        self.assertFalse(utils.names_consistent(mock_names("ABB"), mock_names("ABB")))
        self.assertFalse(utils.names_consistent(mock_names("ABC"), mock_names("ABD")))

def mock_name(name):
    mock = Mock()
    mock.name = name
    return mock

def mock_names(names):
    return [mock_name(name) for name in names]

if __name__ == '__main__':
    unittest.main()