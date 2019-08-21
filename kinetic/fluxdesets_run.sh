#!/bin/bash

ReacList="ACEA, ACEB, ACK, ACN_1, ACN_2, ACS, ATP_MAINTENANCE, ATP_syn, CITRA_SYN, CYTBO, EDA, EDD, ENO, FBA, FBP, FUMA, GDH, GLT, GND, GPM, GROWTH, LPD, MAD, MDH, MQO, NADH_req, NDHII, PCK, PDH, PFK, PGI, PGK, PGL, PIT, PPC, PPS, PTA, PYK, RPE, RPI, SDH, SK, SQR, TPI, XCH_ACE1, XCH_ACE2, XCH_GLC, XCH_P, ZWF"
Field_Separator=$IFS
IFS=,

for reac in $ReacList;
do
	./vmax_run.sh $reac
done

echo "ALL DONE"

