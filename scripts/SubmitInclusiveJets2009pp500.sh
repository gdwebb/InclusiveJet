#!/bin/sh

for runnumber in `cat final_runnumber326.list` ; do
	star-submit-template -template RunInclusiveJetAnalysis.xml -entities runnumber=$runnumber
done
