#!/bin/bash

git clone git@github.com:ijapesigan/fitOU.git
rm -rf "$PWD.git"
mv fitOU/.git "$PWD"
rm -rf fitOU
