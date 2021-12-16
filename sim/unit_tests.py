import unittest
from subprocess import CompletedProcess
from unittest.mock import Mock, mock_open, patch
from sample import Sample
import genome
import sequenced
import processed

import validator

class ValidatorTests(unittest.TestCase):
    # @patch("validator.open", new_callable=mock_open, read_data="mock-text")
    def test_validator(self):
        """ A simple introductory unit test that checks we can run through the
            validator without errors
        """
        validator.shutil.copytree = Mock()
        validator.os.makedirs = Mock()
        validator.os.path.exists = Mock(return_value=False)
        validator.utils.run = Mock()
        validator.glob.glob = Mock(return_value=[''])

        genome.os.path.exists = Mock(return_value=True)

        validator.sequenced.from_results_dir = Mock(return_value=[Mock()])
        validator.processed.from_list = Mock(return_value=[Mock()])
        validator.compare_snps.benchmark = Mock(return_value=[Mock(),{'mock': Mock()}])

        mock_sample = Sample()

        mock_genome = genome.SimulatedGenome('mock-genome', '', '', '', '')
        mock_sample.simulate_genome = Mock(return_value=mock_genome)
        mock_sample.simulate_reads = Mock(return_value=['READ1', 'READ2'])

        validator.performance_test('./', './', [mock_sample])

if __name__ == '__main__':
    unittest.main()