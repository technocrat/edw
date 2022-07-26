#!/bin/bash
cd obj/chunks
cat ../d.csv | parallel --header : --pipe -N2000 'cat >file_{#}.csv'
