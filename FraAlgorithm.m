function output = FraAlgorithm(sig, p)
    output = logssar(sig, p.StimI, p.fs, p.FRA.blankingPeriod);
end