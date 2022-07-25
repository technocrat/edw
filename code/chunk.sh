#!/bin/bash
xsv index d.csv
xsv split splits d.csv
mv splits ../obj
