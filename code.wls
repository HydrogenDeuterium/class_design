#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*pA:=Exp[15.9008-2788.51/(T-52.36)];*)
(*pB:=Exp[16.6076-3295.12/(T-55.6)];*)
(*xA=0.728;xB=1-xA;atm=760;*)
(*root=FindRoot[xA pA+xB pB==atm,{T,273.15}];*)
(*T-273.15/.root*)


(* ::Input:: *)
(*bubblePoint[xA_]:=FindRoot[xA pA+(1-xA) pB==atm,{T,273.15}];*)
(*(T-273.15)/.bubblePoint[0.993]*)
