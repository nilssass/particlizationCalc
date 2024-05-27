(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["xutils`",{"xAct`xTensor`"}]
Unprotect@@Names["xutils`*"];
ClearAll@@Names["xutils`*"];
AddEquation::usage="AddEquation is a book-keeping module, that keeps track of defined equations in the global variable $Equations.
Credits: Alfonso Garc\[IAcute]a-Parrado (http://www.xact.es/xActCourse_Prague/)";
ApplyRule::usage="ApplyRule is used to transform an equation to a rule. 
Credits: Alfonso Garc\[IAcute]a-Parrado (http://www.xact.es/xActCourse_Prague/)";
MoveToLHS::usage="MoveToLHS moves a term located in loc from the right-hand-side of an equation to the left-hand-side.";
MoveToRHS::usage="MoveToRHS moves a term located in loc from the left-hand-side of an equation to the right-hand-side.";
Global`$Equations={};
Begin["`Private`"]
AddEquation[expr_Equal]:=AppendTo[Global`$Equations,Hold[expr]];
Options[ApplyRule]={MetricOn->All,ContractMetrics->True};
ApplyRule[expr_Equal,OptionsPattern[]]:=MakeRule[Evaluate[List@@expr],MetricOn->OptionValue[MetricOn],ContractMetrics->OptionValue[ContractMetrics]]
SetOptions[ApplyRule,MetricOn->All,ContractMetrics->True]
MoveToLHS[eqn_Equal,loc_?IntegerQ]:=Module[{list},
list=List@@eqn;
First@list-(List@@Last@list)[[loc]]==Last@list-(List@@Last@list)[[loc]]]
MoveToRHS[eqn_Equal,loc_?IntegerQ]:=Module[{list},
list=List@@eqn;
First@list-(List@@First@list)[[loc]]==-(List@@First@list)[[loc]]+Last@list]
End[]
Protect@@Names["xutils`*"];
EndPackage[]



