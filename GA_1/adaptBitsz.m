function emaBitsz = adaptBitsz(max_length)

binary = [2,4,8,16,32,64,128,256,512,1028,2056];
emaBitsz = max_length;
for ii = binary
    if emaBitsz < ii
        emaBitsz = find(ii==binary);
        break
    end
end
end