#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:29:26 2024

@author: fabian
"""

import os
import unittest

REPORT = "Rscript {report} {} {} {} {} {} {} {}"

class TestInstallation(unittest.TestCase):
    """
        Test different dependencies to check the installation is correct.
        
        The goal is to idenfy changes in dependencies that may not be correctly
        indicated at creating the conda environment for each
    """
    def test_installation_qc(self):
        import sqanti3_qc
        
    def test_installation_filter(self):
        import sqanti3_filter
        
    def test_installation_rescue(self):
        import sqanti3_rescue
    
    def test_installation_reads(self):
        import sqanti_reads
        
    def test_checkRscript_installation(self):
        assert 0 == os.system("Rscript --version")
    
class TestDependencies(unittest.TestCase):
    
    """
        This batch of test check whether the system-wide dependencies are
        installed and accessible to the user
        
        os.system returns the exit code of the program used. This code is usually 0
        when the program runs correctly, and any positive number when an error
        happened
        
        In particular, shells throw exit code 127 when the command does not 
        exist. Checking if the value is different we can check whether is installed,
        and checking if different from 0 that exited correctly
        
        Checking both values can give more info when an error happened,
        to check if its because it's not installed or if another kind of
        error.        
    """
    
    def test_kallisto(self):
        return_code = os.system("kallisto -h")
        assert 127 != return_code
        assert 0 == return_code
        
    def test_star(self):
        return_code = os.system("star -h")
        assert 127 != return_code
        assert 0 == return_code
    
class TestQC(unittest.TestCase):
    
    def test_fullreport(self):
        command = os.system()
        pass
    
    def test_reportFSM(self):
        pass
    
    def test_reportFSM(self):
        pass
    
    def test_reportFSM(self):
        pass

if __name__ == "__main__":
    unittest.main()