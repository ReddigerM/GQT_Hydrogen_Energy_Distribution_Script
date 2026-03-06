#! /usr/bin/wolframscript

(***************************************************************************************)
(* Article:	On the applicability of Kolmogorov's theory of probability             *)
(*		to the description of quantum phenomena. Part I: foundations           *)
(*                                                                                     *)
(* Journal: 	Quantum Studies: Mathematics and Foundations                           *)
(*                                                                                     *)
(* Author:	Maik Reddiger                                                          *)
(*              Anhalt University of Applied Sciences, Germany                         *)
(* 		maik.reddiger@hs-anhalt.de                                             *)
(*                                                                                     *)
(* DOI:		10.1007/s40509-025-00375-6                                             *)
(*                                                                                     *)
(* Caption:                                                                            *)
(* This wolfram language script computes the numerical data employed for the creation  *) 
(* of Fig. 1 in the article and exports it into the file "gqt_energy_distribution.mx". *)
(* It was written for Mathematica 12.0.0.0 on Linux x86 (64 bit),                      *) 
(* but it may also work for other versions.                                            *)
(* To run it on a machine with Mathematica installed, make the script executable       *) 
(* and execute "./gqt_energy_distribution_script.wls" in a terminal.                   *) 
(* Depending on your machine, executing this script may take a while.                  *) 
(* It will run faster by choosing lower values for nmax, lmax, and mmax.               *)
(***************************************************************************************)

(*******************)
(* ::DEFINITIONS:: *)
(*******************)

(* choice of units: hbar = 1, a0 = 1 and Ry=1 *)

(* set value of Bohr radius to 1 *)
a0 = 1;
(* define hydrogen stationary state basis functions *)
hydrogenpsi[{n_, l_, m_}, {r_, thet_, phi_}] := 
  With[{rho = 2 r/(n a0)}, 
    Sqrt[(2/(n a0))^3 (n - l - 1)!/(2 n (n + l)!)] 
    Exp[-rho/ 2] rho^l LaguerreL[n - l - 1, 2 l + 1, rho] 
    SphericalHarmonicY[l, m, thet, phi]
  ];

(* define energy eigenvalues *)
eeigenv[n_] := -1/n^2;

(* define the trial state via auxiliary functions xi1 and xi2 *)
xi1[r_, thet_, phi_] := hydrogenpsi[{1, 0, 0}, {r, thet, phi}]/Sqrt[2];

xi2[r_, thet_, phi_] := I * hydrogenpsi[{2, 1, 1}, {r, thet, phi}]/Sqrt[2]; 

trialPsi[t_, r_, thet_, phi_] := 
  ( 
  Exp[-I t eeigenv[1]] xi1[r, thet, phi] 
  + Exp[-I t eeigenv[2]]  xi2[r, thet, phi]
  );

(* define frequency (needed for time-evolution) *)
omega = eeigenv[2] - eeigenv[1];

(* define the integrand for the energy random variable *) 
energyintegrand[t_, r_, thet_, phi_] := 
  (
  eeigenv[1] Abs[xi1[r, thet, phi]]^2 
  + (eeigenv[1] + eeigenv[2]) * 
    Re[ Exp[-I t  omega] * Conjugate[xi1[r, thet, phi]] xi2[r, thet, phi]] 
  + eeigenv[2] Abs[xi2[r, thet, phi]]^2
  );


(*******************)
(* ::COMPUTATION:: *)
(*******************)

(* General idea: Compute expectation values for sufficiently many regions in R^3 *)
(*               along with probabilities for those regions.                     *)
(*               Then collect these values to create data, which can be used for *)
(*               making a histogram.                                             *) 

(* :computation of expectation values: *) 

(* number of segments in radial direction (default: 100) *)
nmax = 100;

(* number of segments in polar direction (default: 10) *)
lmax = 10;

(* number of segments in azimuthal direction (default: 20) *)
mmax = 20;

(* this value determines the length of the segments in the radial direction *)
nstretch = 20*nmax/100;

Print["\n Integration intervals in radial direction:"];
Print[Table[{N[Exp[(n - 1)/nstretch] - 1], N[Exp[n/nstretch] - 1]}, {n, nmax}]]; 

Print["\n Integration intervals in polar direction:"];
Print[Table[{Pi (l - 1)/lmax, Pi l/lmax}, {l, lmax}]]; 

Print["\n Integration intervals in azimuthal direction:"];
Print[Table[{2 Pi (m - 1)/mmax, 2 Pi m/mmax}, {m, mmax}]]; 

(* numerical integration *)
probenergy = 
  Table[
    {
      Re[NIntegrate[
        energyintegrand[0, r, thet, phi] *r^2* Sin[thet], 
        {r, Exp[(n - 1)/nstretch] - 1, Exp[n/nstretch] - 1}, 
        {thet, Pi (l - 1)/lmax, Pi l/lmax}, 
        {phi, 2 Pi (m - 1)/mmax, 2 Pi m/mmax}, 
        AccuracyGoal -> 6
      ]], 
      Re[NIntegrate[
        Abs[trialPsi[0, r, thet, phi]]^2 *r^2* Sin[thet], 
        {r, Exp[(n - 1)/nstretch] - 1, Exp[n/nstretch] - 1}, 
        {thet, Pi (l - 1)/lmax, Pi l/lmax}, 
        {phi, 2 Pi (m - 1)/mmax, 2 Pi m/mmax}, 
        AccuracyGoal -> 6
      ]]
    }, 
  {n, nmax}, {l, lmax}, {m,mmax}]; 

(* create list of doubles: *)
(* first entry is the energy expectation value over that region *)
(* second entry is the probability of that region *)
probenergyval = 
  Table[
    {
      Flatten[probenergy, 2][[k, 1]] / Flatten[probenergy, 2][[k, 2]], 
      Flatten[probenergy, 2][[k, 2]]
    }, 
  {k, Length[Flatten[probenergy, 2]]}];

(* this is meant to serve as a check: *) 
(* one should get {-0.625, 1.} *)
Print["\n The total expectation value of the energy 
  and the total probability for the unsorted data is:"]
Print[{Sum[probenergyval[[k, 2]]*probenergyval[[k, 1]], {k, 
   Length[probenergyval]}], 
 Sum[probenergyval[[k, 2]], {k, Length[probenergyval]}]}]

(* compute minimum and maximum values for the energy *)
energymin = Min[Table[probenergyval[[k, 1]], {k, Length[probenergyval]}]]; 
energymax = Max[Table[probenergyval[[k, 1]], {k, Length[probenergyval]}]];

(* compute number of bins (this is a rule of thumb) *)
nbins = Floor[Sqrt[Length[probenergyval]]]; 
Print["\n The number of bins for the data is:"]
Print[nbins]

(* "create empty bins" *)
energybin = 
 Table[
   {
     energymin + (Ceiling[energymax, .0001] - Floor[energymin, .0001])/nbins * (n - 1), 
     0
   }, 
 {n, nbins}]

(* "fill the bins with the data" *)
(
  Do[
    Do[
      If[
        energymin + (Ceiling[energymax, .0001] - Floor[energymin, .0001])/nbins * (n - 1) 
         <= probenergyval[[k, 1]] 
         < energymin + (Ceiling[energymax, .0001] - Floor[energymin, .0001])/nbins * n,
        energybin[[n, 2]] = energybin[[n, 2]] + probenergyval[[k, 2]]
      ], 
    {n, nbins}], 
  {k, Length[probenergyval]}]
)

(* this is meant to serve as a second check: *) 
(* (there may be small deviations from {-0.625, 1.}) *)
Print["\n The total expectation value of the energy 
  and the total probability for the sorted data is:"]
Print[{Sum[energybin[[k, 2]]*energybin[[k, 1]], {k, nbins}], 
 Sum[energybin[[k, 2]], {k, nbins}]}]

(*******************)
(* ::DATA EXPORT:: *)
(*******************)


(* export data in .mx format *)
(* Note: this file format does not use much memory, but it is not cross-platform *)
Export["./gqt_energy_distribution.mx", energybin]; 
Print["Data generation successful."]
