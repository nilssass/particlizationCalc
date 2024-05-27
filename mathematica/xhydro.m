(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["xhydro`",{"xAct`xTensor`"}]
Unprotect@@Names["xhydro`*"];
ClearAll@@Names["xhydro`*"];
CoDiv::usage="CoDiv returns the covariant divergence.";
TFSymmetrize::usage="TFSymmetrize";
Begin["`Private`"]
CoDiv[expr_,tangentM_?VBundleQ,CD_?CovDQ]:=Module[{},(CD[-First@IndicesOf[Free]@expr]@expr)/.First@IndicesOf[Free]@expr->(DummyIn@tangentM)//ScreenDollarIndices]
CoDiv[expr_]:=Module[{tangentM=Global`TangentM,CD=Global`CD},(CD[-First@IndicesOf[Free]@expr]@expr)/.First@IndicesOf[Free]@expr->(DummyIn@tangentM)//ScreenDollarIndices]
TFSymmetrize[expr_,tr_,Delta_,Dim_]:=Module[{},(Symmetrize[expr,IndicesOf[Free]@expr])-(tr) 1/(Dim-1) (Delta[#1,#2]&@@IndicesOf[Free]@expr)];
End[]
Protect@@Names["xhydro`*"];
EndPackage[]



