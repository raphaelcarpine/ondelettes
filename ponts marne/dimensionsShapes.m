if isequal(bridge, 'trilbardou')
    L = 76;
    l = 8.40;
    sensorsPos = [
        14, 8.4
        26, 8.4
        38, 8.4
        14, 0
        26, 0
        38, 0
        ];
    
    shapePlotBridge = @(shape, figTitle) shapePlotPlate([L, l], sensorsPos, shape, [bridge, ' : ', figTitle]);
end