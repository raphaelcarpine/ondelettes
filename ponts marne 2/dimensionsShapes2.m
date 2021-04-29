L = 74 + 2*(0.51+0.125);
l = 8 + 2*0;
sensorsPos = [
    0*L/6, 0
    1*L/6, 0
    2*L/6, 0
    3*L/6, 0
    4*L/6, 0
    5*L/6, 0
    6*L/6, 0
    0*L/6, l
    1*L/6, l
    2*L/6, l
    3*L/6, l
    4*L/6, l
    5*L/6, l
    6*L/6, l
    ];

shapePlotBridge = @(shape, figTitle) shapePlotPlate([L, l], sensorsPos, shape, figTitle);
