def make_table_of_masses(ordered_list_):
    """
       input
           @ordered_list_ : input table produced by function 'get_ordered_list_of_masses()'
       output
           @nice_table : 2d table of mass points that exist
    """
     
    lsp_list = []
    nice_table = []
    nice_table.append('{:^7}'.format(' '))

    for stop in ordered_list_:
        for lsp in ordered_list_[stop]:
            if lsp not in lsp_list:
                lsp_list.append(lsp)
        nice_table.append('{:^7}'.format(stop))

    for lsp in lsp_list:
        nice_table[0] += '{:^7}'.format(str(lsp))

    existence_table = np.zeros((len(ordered_list_), len(lsp_list)), dtype=int)

    for istop, stop in enumerate(ordered_list_):
        for ilsp, lsp in enumerate(lsp_list):
            if lsp in ordered_list_[stop]:
                existence_table[istop][ilsp] += 1

    for irow, row in enumerate(existence_table):
        for entry in row:
            nice_table[irow + 1] += '{:^7}'.format(entry)

    for irow in xrange(len(nice_table)):
        nice_table[irow] += '\n'

    with open('table_of_masses.txt', 'w') as t:
        t.writelines(nice_table)
        t.close()


def get_ordered_list_of_masses(signal_list_):
 
    tree_mass_total_tmp = []
    tree_mass_total = []
    signal_tree_mass_sorted = OrderedDict()

    for signal in signal_list_:
        signal_tree_mass_sorted[signal] = OrderedDict()

        tree_name_mass = [(int(mass.split('_')[1]), int(mass.split('_')[2])) for mass in signal_list_[signal]['trees']]
        tree_mass_total_tmp.append(tree_name_mass)
        tree_name_mass.sort(key=lambda x: int(x[0]))

        for tree in tree_name_mass:
            if str(tree[0]) not in signal_tree_mass_sorted[signal].keys():
                signal_tree_mass_sorted[signal][str(tree[0])] = []

        for tree in tree_name_mass:
            signal_tree_mass_sorted[signal][str(tree[0])].append(tree[1])

        for stop in signal_tree_mass_sorted[signal]:
            signal_tree_mass_sorted[signal][stop].sort(key=lambda x: int(x))        

    for signal in signal_tree_mass_sorted:
        print(signal)

        for stop in signal_tree_mass_sorted[signal]:
            print('   ', stop)

            for lsp in signal_tree_mass_sorted[signal][stop]:
                print('       ', lsp)

    for tree_list in tree_mass_total_tmp:
        for tree in tree_list:
            tree_mass_total.append(tree)

    tree_mass_total.sort(key=lambda x: int(x[0]))
    tree_mass_sorted = OrderedDict()

    for tree in tree_mass_total:
        if tree[0] not in tree_mass_sorted.keys():
            tree_mass_sorted[str(tree[0])] = []

    for tree in tree_mass_total:
        tree_mass_sorted[str(tree[0])].append(tree[1])

    for stop in tree_mass_sorted:
        tree_mass_sorted[stop].sort(key=lambda x: int(x))
       
    return tree_mass_sorted
