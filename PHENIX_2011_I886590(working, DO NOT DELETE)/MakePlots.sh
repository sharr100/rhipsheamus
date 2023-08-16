#!/bin/bash

rivet-mkhtml --pwd --mc-errs Rivet_pp200.yoda
cp -r rivet-plots/ rivet-plots-pp200/

rivet-mkhtml --pwd --mc-errs Rivet_pp62.yoda
cp -r rivet-plots/ rivet-plots-pp62/

