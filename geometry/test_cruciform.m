% test_cruciform
clear classes

pitch       = 2*0.719250;

% Basic surfaces
center.x = pitch/2;
center.y = pitch/2;
angle = 45;
crux = Cruciform(center, angle);

% Basic surfaces
plane_x     = Plane(1, 0, -pitch/2);
plane_y     = Plane(0, 1, -pitch/2);
plane_xyp   = Plane(-1, 1, 0);
plane_xyn   = Plane( 1, 1,-pitch);

% plane_xyn.draw([0 pitch 0 pitch], 100000);
% return
pos_x = SurfaceNode(plane_x, 1);
neg_x = SurfaceNode(plane_x, 0);
pos_y = SurfaceNode(plane_y, 1);
neg_y = SurfaceNode(plane_y, 0);
pos_xyp = SurfaceNode(plane_xyp, 1);
neg_xyp = SurfaceNode(plane_xyp, 0);
pos_xyn = SurfaceNode(plane_xyn, 1);
neg_xyn = SurfaceNode(plane_xyn, 0);
pos_c = SurfaceNode(crux, 1);
neg_c = SurfaceNode(crux, 0);

% Fuel cell
fuelReg1      = Region('fuel1');
fuelReg1.node = OperatorNode(OperatorNode(neg_c, pos_x, Intersection),   pos_xyp, Intersection);
fuelReg2      = Region('fuel2');
fuelReg2.node = OperatorNode(OperatorNode(neg_c, neg_xyp, Intersection), pos_y, Intersection);
fuelReg3      = Region('fuel3');
fuelReg3.node = OperatorNode(OperatorNode(neg_c, neg_y, Intersection),   pos_xyn, Intersection);
fuelReg4      = Region('fuel4');
fuelReg4.node = OperatorNode(OperatorNode(neg_c, neg_xyn, Intersection), pos_x, Intersection);
fuelReg5      = Region('fuel5');
fuelReg5.node = OperatorNode(OperatorNode(neg_c, neg_x, Intersection), neg_xyp, Intersection);
fuelReg6      = Region('fuel6');
fuelReg6.node = OperatorNode(OperatorNode(neg_c, pos_xyp, Intersection), neg_y, Intersection);
fuelReg7      = Region('fuel7');
fuelReg7.node = OperatorNode(OperatorNode(neg_c, pos_y, Intersection), neg_xyn, Intersection);
fuelReg8      = Region('fuel8');
fuelReg8.node = OperatorNode(OperatorNode(neg_c, pos_xyn, Intersection), neg_x, Intersection);
% Moderator
modReg1      = Region('mod1');
modReg1.node = OperatorNode(OperatorNode(pos_c, pos_x, Intersection),   pos_xyp, Intersection);
modReg2      = Region('mod2');
modReg2.node = OperatorNode(OperatorNode(pos_c, neg_xyp, Intersection), pos_y, Intersection);
modReg3      = Region('mod3');
modReg3.node = OperatorNode(OperatorNode(pos_c, neg_y, Intersection),   pos_xyn, Intersection);
modReg4      = Region('mod4');
modReg4.node = OperatorNode(OperatorNode(pos_c, neg_xyn, Intersection), pos_x, Intersection);
modReg5      = Region('mod5');
modReg5.node = OperatorNode(OperatorNode(pos_c, neg_x, Intersection), neg_xyp, Intersection);
modReg6      = Region('mod6');
modReg6.node = OperatorNode(OperatorNode(pos_c, pos_xyp, Intersection), neg_y, Intersection);
modReg7      = Region('mod7');
modReg7.node = OperatorNode(OperatorNode(pos_c, pos_y, Intersection), neg_xyn, Intersection);
modReg8      = Region('mod8');
modReg8.node = OperatorNode(OperatorNode(pos_c, pos_xyn, Intersection), neg_x, Intersection);

g = Geometry(pitch, pitch);
g.addRegion(fuelReg1);
g.addRegion(fuelReg2);
g.addRegion(fuelReg3);
g.addRegion(fuelReg4);
g.addRegion(fuelReg5);
g.addRegion(fuelReg6);
g.addRegion(fuelReg7);
g.addRegion(fuelReg8);

g.addRegion(modReg1);
g.addRegion(modReg2);
g.addRegion(modReg3);
g.addRegion(modReg4);
g.addRegion(modReg5);
g.addRegion(modReg6);
g.addRegion(modReg7);
g.addRegion(modReg8);

g.sketch(4e5);