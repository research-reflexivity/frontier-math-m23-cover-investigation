#!/usr/bin/env python3
"""Negative controls for the calculator-response parser; no network access."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_magma_calculator import parse_response

def response(lines,warning=""):
    return '<calculator><headers><version>2.29-10</version>'+warning+'</headers><results>'+''.join('<line>'+line+'</line>' for line in lines)+'</results></calculator>'

class ResponseTests(unittest.TestCase):
    def test_complete_success(self):
        self.assertEqual(parse_response(response(['PASS_ALL']),['PASS_ALL'])['status'],'PASS')
    def test_wrapped_completion(self):
        self.assertEqual(parse_response(response(['PASS all the','required identities']),['PASS all the required identities'])['status'],'PASS')
    def test_pass_after_assertion_error_is_not_success(self):
        self.assertNotEqual(parse_response(response(['Assertion failed','PASS_ALL']),['PASS_ALL'])['status'],'PASS')
    def test_timeout_even_after_pass_is_not_success(self):
        self.assertNotEqual(parse_response(response(['PASS_ALL'],'<warning>The time limit was exceeded</warning>'),['PASS_ALL'])['status'],'PASS')
    def test_missing_completion_is_not_success(self):
        self.assertNotEqual(parse_response(response(['partial output']),['PASS_ALL'])['status'],'PASS')
    def test_offline_is_not_success(self):
        self.assertNotEqual(parse_response('<calculator><offline>unavailable</offline></calculator>',['PASS_ALL'])['status'],'PASS')
    def test_syntax_error_is_not_success(self):
        self.assertNotEqual(parse_response(response(['Syntax error: bad token','PASS_ALL']),['PASS_ALL'])['status'],'PASS')

if __name__=='__main__':unittest.main()
