import unittest
from subprocess import CompletedProcess
from unittest.mock import Mock, mock_open, patch
from sample import Sample
import genome

import validator

class ValidatorTests(unittest.TestCase):
    def test_validator(self):
        """ A simple introductory unit test that checks we can run through the
            validator without errors
        """
        # Mock I/O on validator
        validator.shutil.copytree = Mock()
        validator.os.makedirs = Mock()
        validator.os.path.exists = Mock(return_value=False)
        validator.utils.run = Mock()
        validator.glob.glob = Mock(return_value=[''])

        validator.compare_snps.benchmark = Mock(return_value=[Mock(),{'mock': Mock()}])

        # Mock Sample
        validator.sequenced.from_results_dir = Mock(return_value=[Mock()])
        validator.processed.from_list = Mock(return_value=[Mock()])
        
        mock_sample = Sample()

        genome.os.path.exists = Mock(return_value=True)
        mock_genome = genome.SimulatedGenome('mock-genome', '', '', '', '')

        mock_sample.simulate_genome = Mock(return_value=mock_genome)
        mock_sample.simulate_reads = Mock(return_value=['READ1', 'READ2'])

        # Run test
        validator.performance_test('./', './', [mock_sample])

if __name__ == '__main__':
    unittest.main()