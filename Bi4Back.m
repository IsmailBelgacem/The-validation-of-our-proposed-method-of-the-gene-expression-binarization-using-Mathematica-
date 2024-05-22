(* ::Package:: *)

(* ::Title:: *)
(*Bi4Back of genetic expression*)


(* ::Author:: *)
(*Ismail Belgacem 2024*)


(* ::Affiliation:: *)
(*Paris*)


BeginPackage["Bi4Back`",{"Booned`"}];


 bin::reg="The variable `1` is neither a positive or a negative regulator of `2`. [ABORT]";


(* ::Input::Initialization:: *)
Bi4Back::usage="Bi4Back[\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"IG_Graph\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"  \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"expressiondata_Association\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"biomarkersphenotypesprofile_Association\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Delta]_Real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Epsilon]_Real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)r\[OHat]le en anglais
PARAMETERS:
\[FilledSquare] role: The binarization of the measured gene expression data. 
\[FilledSquare] Parameters  description
    \[FilledSmallCircle] IG: Interaction Graph
    \[FilledSmallCircle] \!\(\*
StyleBox[\"expressiondata\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"A\",\nFontSlant->\"Italic\"]\)n association of the genes with thiers  measured gene expression data provided.
    \[FilledSmallCircle] \!\(\*
StyleBox[\"biomarkersphenotypesprofile\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"A\",\nFontSlant->\"Italic\"]\)n association all the genes that have a \!\(\*
StyleBox[\"biomarker\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"pheno\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"types\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"profile\",\nFontSlant->\"Italic\"]\).
    \[FilledSmallCircle] \!\(\*
StyleBox[\"\[Delta]_Real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"a\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"value\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\"]\)0.<\[Delta]<=1.\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"used\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"to\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"define\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"the\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"theshould\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"from\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"which\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"a\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"gene\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"  \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"is\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"  \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"considered\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"as\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"expressed\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"or\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"unexpressed\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)
    \[FilledSmallCircle] \!\(\*
StyleBox[\"\[Epsilon]_Real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"a\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"real\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"value\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\"]\)0.<\[Epsilon]<=1.)used to define \!\(\*
StyleBox[\"the\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"theshould\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"from\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"which\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)regulators are considered to have almost a very similar genes expressions.
"


(* ::Text:: *)
(*\[FilledSquare] Inconsistent Regulators function: *)
(*\[FilledSmallCircle] Definition: Here, the concept  of inconsistency  is defined according to the current state of the considered  target.  The inconsistent regulators are all the regulators that could not be in any cases  responsible for providing  the binarized value of the  target. However, if all the regulators of the considered target are inconsistent, then there exists a confusion according the interaction graph and in this case  the inconsistency becomes stronger. Where, at least one regulator should be consistent  and which is responsible for the current Boolean state  of the target. To avoid therefore  this problem,  all the regulators of the considered target are substituted by not assigned values. *)
(*\[FilledSquare] The function role: Select the already binarized  inconsistent regulators of the considered target. *)
(*\[FilledSquare] Parameters  description*)
(*    \[FilledSmallCircle] target: is the considered target.*)
(*    \[FilledSmallCircle]  ig: Interaction Graph*)
(*    \[FilledSmallCircle] binvar_Association: An association of all the  genes that have already binarized values (binarized according to the provided gene expression data,  biomarkers phenotypes profiles  or with a forward and back propagation). *)


InconsistentRegulators[target_,ig_Graph, binvar_Association]:=Select[Regulators[ig,target,{1}], Function[reg, (KeyExistsQ[binvar,reg] \[And] ((binvar[reg] \[And] \[Not] binvar[target]) || (\[Not]binvar[reg] \[And] binvar[target])))]]\[Union]Select[Regulators[ig,target,{-1}], Function[reg, (KeyExistsQ[binvar,reg] \[And] ((\[Not]binvar[reg] \[And] \[Not] binvar[target]) || (binvar[reg] \[And] binvar[target])))]]


(* ::Text:: *)
(*\[FilledSquare]  GetValueBackPropagation function : *)
(*\[FilledSquare] role:  A back  propagation  of the Boolean value  of the  considered already binarized target to  provide the Boolean value of the regulator  that is picked out  or defined (from all the regulators of the target)  as the one that has  the most expressed value (if it is an activator ) or the most unexpressed value  (if is an inhibitor). This means define  the Boolean value of the regulator that is selected according the expression data as the one which is  responsible for the currant value of  target.*)
(*\[FilledSquare] Parameters  description*)
(*    \[FilledSmallCircle] t:  The regulator that is defined (according the expression data) as  the one which  is  responsible for the current value of  target.*)
(*    \[FilledSmallCircle] ig: Interaction Graph*)
(*    \[FilledSmallCircle] target: is the considered target.*)
(*    \[FilledSmallCircle] boolvalofg: is the Boolean value of the considered target.*)


GetValueBackPropagation[t_, target_,ig_Graph,boolvalofg:(True|False)]:=
Which[
MemberQ[Regulators[ig,target,{1}],t],   boolvalofg,
MemberQ[Regulators[ig,target,{-1}],t],\[Not] boolvalofg,
True, Message[bin::reg, t,target];Abort[]
];


(* ::Text:: *)
(*\[FilledSquare]  Harmonisation function : *)
(*\[FilledSquare] role: Binarize another regulator of the considered target that  has almost a similar  gene expression value compared to  the regulator that is already picked out as the one which is  responsible for the current value of the target.*)
(*\[FilledSquare] Parameters  description*)
(*    \[FilledSmallCircle] t: Another regulator of the considered target that has almost the same  gene expression amount compared to  the regulator already picked out as the one which is  responsible for the current value of  target.*)
(*    \[FilledSmallCircle] z: The regulator already selected as the one which is  responsible for the current value of  target (the binarized value of output of function before GetValueBackPropagation).*)
(*    \[FilledSmallCircle] ig: Interaction Graph*)
(*    \[FilledSmallCircle] target: The considered target.*)
(*    \[FilledSmallCircle] boolvalofz: The Boolean value of the regulator already selected as the one which is  responsible for providing  the current Boolean value of the  target.*)
(*    *)


(* ::Input::Initialization:: *)
Harmonisation[t_,z_, target_,ig_Graph,boolvalofz:(True|False)]:=
Which[
(MemberQ[Regulators[ig,target,{1}],t]&& MemberQ[Regulators[ig,target,{1}],z])||(MemberQ[Regulators[ig,target,{-1}],t]&& MemberQ[Regulators[ig,target,{-1}],z]), boolvalofz,
(MemberQ[Regulators[ig,target,{-1}],t]&& MemberQ[Regulators[ig,target,{1}],z])||(MemberQ[Regulators[ig,target,{1}],t]&& MemberQ[Regulators[ig,target,{-1}],z]),\[Not]boolvalofz,
True, Message[bin::reg, t,target];Abort[]
];


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Main Function*)


(*Begin["`Private`"]*)


(* ::Input::Initialization:: *)
Bi4Back[ ig_Graph,  Expressiondata_Association, biomarkersphenotypesprofile_Association,\[Delta]_Real/; 0.<\[Delta]<=1.,\[Epsilon]_Real/; 0.<\[Epsilon]<=1.]:= Module[{positiveregs,negativeregs,R,\[Tau]values,chosenregs,undefinedvars, previousbinarizedvar,alreadybinarizedvars, binarizedvar,existsupdate,  genename,int,expressiondata,Cm},
     binarizedvar=<||>;
genename=Keys[Expressiondata]\[Union] Complement[Keys[biomarkersphenotypesprofile],Keys[Expressiondata]];
int= DeleteCases[Expressiondata,"Na"|"Nan"|"Missing"|"Nothing"|""|"None"|Na|Nan|Missing|Nothing|None];
If[Max[int]>=1,int=N[(int-Min[int])/(Max[int]-Min[int]) ]];(*Echo[int];*)
Cm=Complement[VertexList[ig],Keys[int]];
If[Cm =!={}, expressiondata=Join[int, AssociationThread[Cm,0.5]]
,expressiondata=int];
(* Initialization of the binarisation profile *)
Do[
Which[
(KeyExistsQ[biomarkersphenotypesprofile,var])\[And](\[Not]MissingQ[biomarkersphenotypesprofile[var]]), binarizedvar[var]=biomarkersphenotypesprofile[var],
expressiondata[var]<\[Delta], binarizedvar[var]=False,
expressiondata[var]>=1-\[Delta], binarizedvar[var]=True];
,{var,genename}];

existsupdate=True;
i=1;

While[existsupdate \[And] i<=(Length[genename]+1) ,
 existsupdate=False;
	(* Forward propagation *)
previousbinarizedvar=(genename/.binarizedvar);
undefinedvars=Complement [genename, Keys[binarizedvar]];
Do[
positiveregs =Regulators[ig,var,{1}];
negativeregs =Regulators[ig,var,{-1}];
Which[
       (* all the activators are true  and all the inhibitors are false *) 
Regulators[ig,var]=!={}\[And] AllTrue[positiveregs/.binarizedvar,#===True &] \[And]  AllTrue[negativeregs/.binarizedvar, #===False &], 
binarizedvar[var]=True,

        (* all the activators are false  and all the inhibitors are true *) 
Regulators[ig,var]=!={}\[And]  AllTrue[positiveregs/.binarizedvar,#===False &] \[And]  AllTrue[negativeregs/.binarizedvar, #===True &],
binarizedvar[var]=False
]
,{var,undefinedvars}];

(* Backward propagation *)

alreadybinarizedvars=Keys[binarizedvar];
Do[
If[(KeyExistsQ[binarizedvar,var]),
R=Complement[Regulators[ig,var],InconsistentRegulators[var,ig,binarizedvar]];
Which[
R==={} \[And] Regulators[ig,var]=!={}, KeyDropFrom[binarizedvar,Flatten[{Regulators[ig,var],var}]],
R=!= {} \[And] Regulators[ig,var]=!={} \[And] AllTrue[R,\[Not]KeyExistsQ[binarizedvar,#]&],
\[Tau]values=<||>;
Do[
Which[
MemberQ[Regulators[ig,var,{1}],reg],\[Tau]values[reg]=Boole[binarizedvar[var]](1-expressiondata[reg])+ expressiondata[reg](1-Boole[binarizedvar[var]]) ,MemberQ[Regulators[ig,var,{-1}],reg],\[Tau]values[reg]=Boole[binarizedvar[var]]*expressiondata[reg] + (1-expressiondata[reg])*(1-Boole[binarizedvar[var]])
,True, \[Tau]values[reg]=Infinity
];
,{reg,R}];

chosenregs=Keys@MinimalBy[\[Tau]values,Identity];
Do[
binarizedvar[creg]=GetValueBackPropagation[creg, var,ig, binarizedvar[var]];AppendTo[alreadybinarizedvars,creg]
,{creg,chosenregs}];

Do[
If[ 
Abs[\[Tau]values[reg]-\[Tau]values[First[chosenregs]]]<= \[Epsilon], binarizedvar[reg]=Harmonisation[reg,First[chosenregs], var,ig, binarizedvar[First[chosenregs]]];AppendTo[alreadybinarizedvars,reg]
]
,{reg,Complement[R,chosenregs]}];
];
]
,{var,alreadybinarizedvars}];
existsupdate= (genename/.binarizedvar)=!=previousbinarizedvar;
i=i+1;(*Echo[i]; Echo[binarizedvar];*)
] ;(* end global while *)
BooleanProfile = binarizedvar
] 


(* ::Input::Initialization:: *)
(*End[];*)


(* ::Input::Initialization:: *)
EndPackage[]
