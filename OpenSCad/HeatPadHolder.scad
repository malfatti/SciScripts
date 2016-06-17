PadLen = 15.3;
PadWid = 6;
PadHei = 1;

union(){
    cube([PadLen, PadWid, PadHei]);
    translate([1, 0, 1])cube([1,1,1]);
}