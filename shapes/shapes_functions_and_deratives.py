def shape_functions_create(ksi, eta):
    N = []
    temp_count = (1 / 4) * (1 - ksi) * (1 - eta)
    N.append(temp_count)
    temp_count = (1 / 4) * (1 + ksi) * (1 - eta)
    N.append(temp_count)
    temp_count = (1 / 4) * (1 + ksi) * (1 + eta)
    N.append(temp_count)
    temp_count = (1 / 4) * (1 - ksi) * (1 + eta)
    N.append(temp_count)
    return N


def dn_dksi(eta):
    dn_ksi = []
    temp_count = (-1 / 4) * (1 - eta)
    dn_ksi.append(temp_count)
    temp_count = (1 / 4) * (1 - eta)
    dn_ksi.append(temp_count)
    temp_count = (1 / 4) * (1 + eta)
    dn_ksi.append(temp_count)
    temp_count = (-1 / 4) * (1 + eta)
    dn_ksi.append(temp_count)
    return dn_ksi


def dn_deta(ksi):
    dn_eta = []
    temp_count = (-1 / 4) * (1 - ksi)
    dn_eta.append(temp_count)
    temp_count = (-1 / 4) * (1 + ksi)
    dn_eta.append(temp_count)
    temp_count = (1 / 4) * (1 + ksi)
    dn_eta.append(temp_count)
    temp_count = (1 / 4) * (1 - ksi)
    dn_eta.append(temp_count)
    return dn_eta
