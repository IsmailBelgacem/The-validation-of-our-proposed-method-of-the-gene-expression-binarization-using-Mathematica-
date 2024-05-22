(* ::Package:: *)

(* ::Title:: *)
(*Other Functions needed :*)


(* ::Author:: *)
(*Ismail Belgacem 2024*)


(* ::Affiliation:: *)
(*Paris*)


BeginPackage["Functionsneededpackage`",{"Booned`"},{"Bi4Back`" }];


(* ::Section:: *)
(*Hills function*)


n=60


Clear[hillp]
hillp[x_Real,thr_Real]:=  x^n/(thr^n+x^n);


Clear[hillm]
hillm[x_Real,thr_Real]:=  thr^n/(x^n+thr^n);


(* ::Section:: *)
(* Random thresholds function*)


Clear[Thresholding]
Thresholding[vars_List,basicthreshold_, dispersion_]:=AssociationThread[ vars, Array[basicthreshold+RandomReal[{-dispersion,dispersion}]&, Length[vars]]]


(* ::Section:: *)
(*A Function that extract three experiments  from an  ODE system simulation.*)


Clear[ThreeExperienceFromODE]
ThreeExperienceFromODE::usage="ThreeExperienceFromODE[F,system, tf, StableState, \[Theta],\[CapitalDelta]] extract three experiments from  an  ODE system simulation. 
PARAMETERS:
\[FilledSmallSquare] F : Boolean network.
\[FilledSmallSquare] system : A system of ODE equations.
\[FilledSmallSquare] tf : The final time of simulation.
\[FilledSmallSquare] StableState: The stady state that correspends  to the initial condition inclouded in the ODE system. 
\[FilledSmallSquare] \[Theta] : Threshold of expression (Association var \[Rule] Real)
\[FilledSmallSquare] \[CapitalDelta]: precision of values (default \!\(\*SuperscriptBox[\(10\), \(-10\)]\))
"
ThreeExperienceFromODE[F_ ,system_List, tf_, StableState_Association, \[Theta]_Association]:=
Module[{Tm,Sout,Sol, Reg,M,bin,idstable,Vec,i,sol, firststable, Exp, TmExp},
Tm=Table[t,{t,0,tf,0.01}];
(* system resolution*) 
Sout=NDSolve[system,Keys[F],{t,0,tf}]; 

(* time step discretization (time step 0.01) *)
Sol=Chop[SetPrecision[Table[First[Evaluate[Map[#[t]&, Keys[F]]/. Sout]],{t,0,tf,0.01}],3],10^-4];

(* Boolean profile for any step w.r.t. the threshold *)
Reg=Table[Association@Table[M=AssociationThread[Keys[F],sol] ; var->M[var]> \[Theta][var],{var,Keys[F]}], {sol,Sol}];

(* Find the stable region *) 
idstable=Table[  If[AllTrue[Keys[FB],(Reg[[i]][#]==StableState[#]&)],i,Nothing],{i,1,Length[Reg]} ];
If[idstable=!={},
(* Take 3 positions in the stable region *)
If[Length[idstable]>3,
Vec={idstable[[1]]+2, idstable[[ Round[Length[idstable]/2]]]  ,idstable[[-1]] }, Vec=idstable];

(* OR Vec={l[[1]]+1,l[[1]]+10, l[[1]]+Floor[SetPrecision[Divide[(l[[Length[l]]]-l[[1]]),2],3]] }*)

Exp=Table[AssociationThread[Keys[F],Sol[[i]]],{i, Vec}];
TmExp=AssociationThread[Tm[[Vec]],Exp], TmExp={}];
TmExp
(* MinMax Normalization  *)
(*Exp2=Exp1;
Do[If[Max[Exp1[[i]]]>=1,Exp1[[i]]=(Exp1[[i]]-Min[Exp1[[i]]])/(Max[Exp1[[i]]]-Min[Exp1[[i]]]) ], {i,Length[Exp1]}];
TmExp1=AssociationThread[Tm[[Vec]],Exp1]; Exp2=AssociationThread[Tm[[Vec]],Exp2],TmExp1={}];
{Exp2,TmExp1}*)
]



(* ::Section:: *)
(*A function that binaries the three experiments extracted.*)


Clear[ThreeExperimentsBinarisation];


ThreeExperimentsBinarisation::usage="ThreeExperimentsBinarisation[F_, Ex_List , StableState_] binarizing  the three experiments and determining the distance from the stable state. 
PARAMETERS:
\[FilledSmallSquare] F : Interaction Graph.
\[FilledSmallSquare] Ex : The  experiments values  (Association var \[Rule] Real)
\[FilledSmallSquare] StableState: The corresponding stady state. 
"


Clear[ThreeExperimentsBinarisation]
ThreeExperimentsBinarisation[G_Graph, Ex_List , StableState_]:=Module[{Biomarker,genenames,Binarizedgenes,BolProfiles1,input,Bi,Tab1,HTH,Bin,BinR,BolProfile2,St,Dis},
Biomarker=<||>;
genenames=VertexList[G];
Binarizedgenes={};
BolProfiles1={};
Do[
Bi=Bi4Back[G,input,Biomarker, 0.15, 0.05];
(*Echo[Bi];*)
AppendTo[Binarizedgenes,Bi];
Tab1=AssociationThread[genenames->None];
HTH=Append[Bi,KeyComplement[{Tab1,Bi}]];
Bin=KeySortBy[HTH,Position[genenames,#]&];
AppendTo[BolProfiles1, Bin];
,{input, Ex}];
BolProfiles1
]


(* ::Section:: *)
(*A function that  calculates a distance from the stable state.*)
(*(A distance from the stable state is determined after binarizing the experiments)*)


(* ::Input:: *)
(*Distanceofdissimilarity::usage="Distanceofdissimilarity[Res_List , StableState_] determines the distance from the stable state. *)
(*PARAMETERS:*)
(*\[FilledSmallSquare] Res : The  experiments biniry profiles  *)
(*\[FilledSmallSquare] StableState: The corresponding stady state. *)
(*"*)


Clear[Distanceofdissimilarity]
Distanceofdissimilarity[ Res_List , StableState_]:=
Module[{d},
Table[If[Values[var]===(Keys[var]/.StableState), d=0, d=0; Do[If[(in/.var)=!=(in/.StableState),d=d+1], {in, Keys[var]}]; d=(d/Length[var])];d,{var,Res}]
]


(* ::Section:: *)
(*A function that  display a nice graph of interaction graph.*)


Clear[IGraphView]
IGraphView[igf_Graph,highlights_:{}]:=DynamicModule[{All1,arrowsize=0.015, graphview,valcentrality,nodecentrality,location={Left,Bottom},style=Directive[FontSize->11, FontFamily->"Yu Gothic UI Semibold"]},
All1[g_]:=ConstantArray[1,VertexCount[g]];
Highlighted[Manipulate[
valcentrality=funcentrality[igf];
nodecentrality=Association@Thread[ VertexList[igf]-> valcentrality/Total[valcentrality]];
Highlighted[
graphview=HighlightGraph[SetProperty[If[self,igf,EdgeDelete[igf,EdgeList[igf, $v_ \[DirectedEdge] $v_]]]
, { 
GraphLayout->layout,
ImageSize->imagesize, 
VertexSize->Normal[ vscale (0.025+nodecentrality)],
VertexLabelStyle-> style,
EdgeShapeFunction-> GraphElementData[{"Arrow","ArrowSize"->arrowsize}]}]
,First/@highlights
,Epilog-> If[highlights=={},{},Inset[SwatchLegend[Last/@highlights[[;;,1]], highlights[[;;,2]], LabelStyle->style] ,location,location]]
]
,Frame->True, Background->GrayLevel[1], FrameStyle->Directive[ GrayLevel[0.3]]]
,{{self,False,"Self"},{True,False}}
,{{layout,Automatic,""},{Automatic,"LayeredDigraphEmbedding","SpringElectricalEmbedding","SpringEmbedding","CircularEmbedding","RadialEmbedding"}}
, {{funcentrality,All1,"Centrality"},{All1  -> "None", DegreeCentrality-> "Degree",  BetweennessCentrality-> "Betweenness",  ClosenessCentrality-> "Closeseness", EigenvectorCentrality-> "Eigenvector"}, ControlType->PopupMenu}
,{{vscale,3,"\[FilledCircle] Size"},0,10,0.5, ControlType->SetterBar}
,{{imagesize,1000,"Img. Size"},{Large,1000,1300,Full}}
, ControlPlacement->{Top,Top,Top,Top,Top,Bottom,Bottom}
,LabelStyle->style
,ContinuousAction->False
, Paneled->False
], Frame->True, Background->GrayLevel[0.84], FrameStyle->{ GrayLevel[0.3]}, RoundingRadius->5]]


EndPackage[]
