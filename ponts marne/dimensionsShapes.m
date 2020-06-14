if isequal(bridge, 'trilbardou')
    L = 76;
    l = 8.40;
    sensorsPos = [13, 0; 13, 8.4; 26, 0; 26, 8.4; 38, 0; 38, 8.4];
    
    shapePlotBridge = @(shape, figTitle) shapePlotPlate([L, l], sensorsPos, shape, [bridge, ' : ', figTitle]);
end