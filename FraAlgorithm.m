function output = FraAlgorithm(sig, p)
    output = logssar(sig, p.StimI, p.fs, 1e-3);
end