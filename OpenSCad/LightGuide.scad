Thick = 4; CircleRes = 500;
CupDiam = 7.9; CupLen = 13.5;
HolderSide = 42.5; ScrewDiam = 2.5;

OutHolderSide = HolderSide + Thick;
HolderCenter = (OutHolderSide/2) - (HolderSide/2);

OutCupDiam = CupDiam + Thick*2;
OutCupCenter = OutHolderSide/2;
CupCenter = HolderSide/2;

HoleSide = HolderSide/2 - Thick;
Hole1Center = Thick;
Hole2Center = HolderSide - HoleSide;

Screw1X = OutHolderSide/2; Screw1Z = Thick;// + ScrewDiam/2;
Screw2Y = CupCenter-CupDiam/2; Screw2Z = CupLen;

module Holder(){
    module Outer(){
        union(){
            cube([HolderSide+Thick, HolderSide+Thick, Thick+Thick]);
            translate([OutCupCenter, OutCupCenter, Thick]) cylinder(r=OutCupDiam/2, h=CupLen, $fn=CircleRes);
        }
    }

    module Inner(){
        union(){
            translate([HolderCenter, HolderCenter, 0]) cube([HolderSide, HolderSide, Thick]);
            translate([OutCupCenter, OutCupCenter, Thick]) cylinder(r=CupDiam/2, h=CupLen, $fn=CircleRes);
        }
    }
    
    difference(){Outer(); Inner();}
}

module Holes(){
    translate([Hole1Center, Hole1Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole2Center, Hole2Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole1Center, Hole2Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole2Center, Hole1Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    
    translate([Screw1X, Thick/2, Screw1Z]) rotate([90,0,0]) cylinder(r=ScrewDiam/2, h=Thick/2, $fn=CircleRes);
    translate([Screw1X, Screw2Y, Screw2Z]) rotate([90,0,0]) cylinder(r=ScrewDiam/2, h=Thick, $fn=CircleRes);
    
}
difference(){Holder(); Holes();}