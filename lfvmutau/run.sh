#!/bin/bash

# Analyze the H->mutau channel

set -o nounset
set -o errexit

export jobid=SubmitMuTauSingleMU
rake mttight


