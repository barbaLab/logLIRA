function output = FraAlgorithm(sig, p)
    output = logssar(sig, p.StimI, p.fs, 1e-3, 'SaturationVoltage', 6);
end