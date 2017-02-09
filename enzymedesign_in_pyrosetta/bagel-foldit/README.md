# Foldit model of BglB for the Bagel project

This repository contains a macromolecular model of BglB based on the crystal structure 2JIE that is used for design efforts in the Siegel group's Bagel project. Read on to learn the ins and outs of using Foldit for enzyme design.

## Load your protein 

Start ---> Intro Puzzles --> Level 1-1 -->  Intro Puzzle loaded hit:  Ctrl-Alt-Shift-A

Select all four files to load into Foldit. 

## Use the advanced design interface

Menu --> Selection Interface
Menu --> General Options --> Show Advanced GUI
Menu --> View Options: 

## View options 

+ "Cartoon Thin"
+ "Score/Hydro+CPK"
+ show/hide bonds/voids as needed for context 

## Save Your Protein 

Ctrl-Alt-Shift-S

## Residue energy

Mouse over residue so it is highlighted and hit Tab. 

This info panel gives a breakdown of the various score terms for the selected residue. Importance of score is in order as it appears in the menu

Viewing operations
------------------

+ Zoom to Area of Interest: Mouse over residue of interest so it is highlighted;  hit Shift-Q  (no click) I often move around until I see the ligand and then do this

+ BackShading:  Ctrl-Shift-Click (on black backdrop) and drag

+ Front Clipping:  Ctrl-Alt-Click (on black backdrop) and drag

+ Measuring:  Mouse over atom1 of interest, hit Ctrl-Alt-Shift-C (no click), repeat for atom 2 (reports distance), atom 3 (angle), and atom 4 (dihedral)


Selection operations 
--------------------

+ See upper left corner of Foldit

+ SphereAround:  Ctrl-Shift-Click on residue and drag. 


Using the Lua Terminal  
----------------------

+ Main --> Script Terminal
+ Select specific residue:  selection.Select(<res#>)
+ Zoom to specific residue:  ui.CenterViewport(<res#>)

Design operations
-----------------
+ Click residue and hit mutate (m).  Selecting a subset of residues (lower right boxes) can be used to do directed designs.  This selection will be used in both mutate AND shake (i.e. shake now mutates those residues while everything else stays native).


Rules of Thumb
--------------

+ **Never** wiggle (w) without selecting a sub-selection!
+ Sometimes high energies (mutations/conformations) are need to get to the low energy
+ Native amino acid or close to native amino acid is generally better
+ Manually force the interaction you want and then side-chain wiggle (E) those residues and the ones around it
+ General order of operations:  E --> S --> E
+ If score is not dropping faster than 0.05/s you are likely at an energy minimum. This is true **only for** wiggling (w/e).  If you are shaking (s) pay attention to the number of cycles on the upper left corner.  Often the low energy is found after 5-10 cycles (cycles is the number in the parenthesis after the timer).

General Energy Rules (on a per amino acid basis)
------------------------------------------------
+ Energies less than 0 are good.
+ Energies 0 - +2 can be ok if clashing (fa_rep) is low (i.e. <1)
+ Energies +2 - +5 are iffy… but can work if you think the "static" picture is missing the enabling feature
+ Energies +5 - +10 are are unlikely… something will likely have to move for this to work
+ Energies >+10 are are really unlikely to work… but go for it if you have a good reason
