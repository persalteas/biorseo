import json
import numpy as np
import matplotlib.pyplot as plt
import os.path

# Creates a violin plot of the distribution of the number of 'motifs' in the Isaure pattern library
def stats_library():
    with open('/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/motifs_28-05-2021.json') as f:
        data = json.load(f)

        nb_motifs = json.dumps(data).count("sequence")
        print("nombre de motifs: " + str(nb_motifs))

        tab_seq_length = []
        tab_seq_length_with_delimiter = []
        for i in range(1007):
            test = str(i) in data
            if test:
                sequence = data[str(i)]["sequence"]
                count_delimiter = sequence.count('&')
                tab_seq_length.append(len(sequence) - count_delimiter)
                tab_seq_length_with_delimiter.append(len(sequence))

        data_to_plot = [np.array(tab_seq_length), np.array(tab_seq_length_with_delimiter)]
        min1 = np.amin(data_to_plot[0])
        max1 = np.amax(data_to_plot[0])
        median1 = np.median(data_to_plot[0])

        min2 = np.amin(data_to_plot[1])
        max2 = np.amax(data_to_plot[1])
        median2 = np.median(data_to_plot[1])
        fig = plt.figure()

        ax = fig.add_axes([0, 0, 1, 1])

        label1 = "nombre de nucléotides" + "\n minimum: " + str(min1) + "    mediane: " + str(
            median1) + "    maximum: " + str(max1)
        label2 = "nombre de nucléotides + nombre de &" + "\n minimum: " + str(min2) + "    mediane: " + str(
            median2) + "    maximum: " + str(max2)
        labels = [label1, label2]
        ax.set_xticks(np.arange(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        ax.set_ylabel('longueurs des motifs')
        ax.set_xlabel('motifs')

        violins = ax.violinplot(data_to_plot, showmedians=True)

        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vp = violins[partname]
            vp.set_edgecolor('black')
            vp.set_linewidth(1)

        for v in violins['bodies']:
            v.set_facecolor('red')
        plt.title("Répartition des longueurs des motifs dans la bibliothèque (" + str(nb_motifs) + " motifs)")
        plt.savefig("statistiques_motifs_Isaure.png", bbox_inches='tight')

# Returns the list in half
def get_half(list_name):
    first_half = []
    second_half = []
    if len(list_name) % 2 == 0:
        middle = len(list_name) / 2
    else:
        middle = len(list_name) / 2 + 0.5

    for i in range(int(middle)):
        first_half.append(list_name[i])

    for i in range(int(middle)):
        if i + int(middle) < len(list_name):
            second_half.append(list_name[i + int(middle)])

    return [first_half, second_half]

# Returns the list of name of the sequence in the benchmark.dbn
def get_list_name_bm():
    path_file = '/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/benchmark.dbn'
    my_file = open(path_file, "r")
    list_name = []
    count = 0
    for line in my_file:
        if count % 4 == 0:
            list_name.append(line[5:len(line) - 2])
        count = count + 1
    my_file.close()

    return list_name

# Returns a 2d array containing for each sequence of the benchmark the number of 'motifs' inserted of each solution
def get_nb_motifs_by_seq(type_file):
    list_name = get_list_name_bm()
    tab = []
    for name in list_name:
        path_file = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/results/test_" + name + type_file
        if os.path.exists(path_file):
            list_nb_motifs = []
            my_file = open(path_file, "r")
            name = my_file.readline()
            seq = my_file.readline()
            count = 0
            for line in my_file:
                if count % 2 == 0:
                    tab_split = line.split("+")
                    list_nb_motifs.append(len(tab_split) - 1)
                count = count + 1
            my_file.close()
            tab.append(list_nb_motifs)
    return tab

# Creates a violin plot that shows the distribution of the number of patterns per solution for each sequence of the benchmark
def stats_nb_motifs_in_result(type_file):
    list_name = get_list_name_bm()
    tab = get_nb_motifs_by_seq(type_file)

    list_median_str = []
    for i in range(len(tab)):
        list_median_str.append(np.median(tab[i]))

    tab = [x for _, x in sorted(zip(list_median_str, tab))]
    list_name = [x for _, x in sorted(zip(list_median_str, list_name))]

    if (len(tab) % 2 == 0):
        absciss = len(tab) / 2
    else:
        absciss = len(tab) / 2 + 0.5

    divide_name = get_half(list_name)
    divide_tab = get_half(tab)

    arr = np.array(tab)
    fig, ax = plt.subplots()
    plt.figure(figsize=(15, 4), dpi=200)
    plt.xticks(rotation=90)
    violins = plt.violinplot(divide_tab[0], showmedians=True)
    for i in range(int(len(divide_tab[0]))):
        y = divide_tab[0][i]
        x = np.random.normal(1 + i, 0.04, size=len(y))
        plt.scatter(x, y)
        plt.xticks(np.arange(1, len(divide_tab[0]) + 1), divide_name[0])

    plt.xlabel('nom de la séquence')
    plt.ylabel('nombre de motifs insérés dans la structure prédite')
    plt.title("Répartition du nombre de motifs insérés par résultat pour chaque séquence")

    for part in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[part]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor('black')

    plt.savefig("statistiques_nb_motifs_inseres_Isaure_" + type_file[8:] + "_1.png", bbox_inches='tight')

    plt.figure(figsize=(15, 4), dpi=200)
    plt.xticks(rotation=90)
    violins = plt.violinplot(divide_tab[1], showmedians=True)
    for i in range(int(len(divide_tab[1]))):
        y = divide_tab[1][i]
        x = np.random.normal(1 + i, 0.04, size=len(y))
        plt.scatter(x, y)
        plt.xticks(np.arange(1, len(divide_tab[1]) + 1), divide_name[1])

    plt.xlabel('nom de la séquence')
    plt.ylabel('nombre de motifs insérés dans la structure prédite')
    plt.title("Répartition du nombre de motifs insérés par résultat pour chaque séquence (2ème partie)")

    for part in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[part]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor('black')
    plt.savefig("statistiques_nb_motifs_inseres_Isaure_" + type_file[8:] + "_2.png", bbox_inches='tight')

# Returns the grouping of the number of inserted 'motif' for all solutions of all sequences of the benchmark
# according to the extension in argument
def get_all_nb_motifs_by_type(type_file):
    tab = get_nb_motifs_by_seq(type_file)
    tab_all = []
    for i in range(len(tab)):
        for j in range(len(tab[i])):
            tab_all.append(tab[i][j])
    return tab_all

# Create a figure containing the violin plot for MEA + function E, MEA + function F, MFE + function E and MFE + function F
# Each violin plot show the distribution of the number of inserted 'motif' by solution
def stats_nb_motifs_all():
    list_name = get_list_name_bm()
    tab_all_E_MEA = get_all_nb_motifs_by_type(".json_pmE_MEA")
    tab_all_F_MEA = get_all_nb_motifs_by_type(".json_pmF_MEA")
    tab_all_E_MFE = get_all_nb_motifs_by_type(".json_pmE_MFE")
    tab_all_F_MFE = get_all_nb_motifs_by_type(".json_pmF_MFE")

    data_to_plot = [tab_all_E_MEA, tab_all_F_MEA, tab_all_E_MFE, tab_all_F_MFE]
    fig = plt.figure()
    fig.set_size_inches(6, 3)
    ax = fig.add_axes([0, 0, 1, 1])

    labels = ['MEA + E', 'MEA + F', 'MFE + E', 'MFE + F']
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Répartition du nombre de motifs insérés par résultat pour chaque séquence")
    ax.set_ylabel("nombre de motifs insérés dans la structure prédite")

    violins = ax.violinplot(data_to_plot, showmedians=True)

    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[partname]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor("black")
    plt.savefig('repartition_nb_motifs.png', dpi=200, bbox_inches='tight')

# Returns the list of the number of solutions and the list of the length of each sequence of the benchmark
def get_nb_solutions_and_sizes_by_seq(type_file):
    list_name = get_list_name_bm()
    list_nb_solutions = []
    list_size = []
    for name in list_name:
        path_file = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/results/test_" + name + type_file
        if os.path.exists(path_file):
            my_file = open(path_file, "r")
            name = my_file.readline()
            seq = my_file.readline()
            list_size.append(len(seq))
            count = 0
            nb = 0
            for line in my_file:
                if count % 2 == 0:
                    nb = nb + 1
                count = count + 1
            list_nb_solutions.append(nb)
            my_file.close()
    return [list_nb_solutions, list_size]

# Creates 4 violin plots that shows the distribution of the number of solutions for each sequence of the benchmark
def stats_nb_solutions():
    list_name = get_list_name_bm()
    tab_all_E_MEA = get_nb_solutions_and_sizes_by_seq(".json_pmE_MEA")[0]
    tab_all_F_MEA = get_nb_solutions_and_sizes_by_seq(".json_pmF_MEA")[0]
    tab_all_E_MFE = get_nb_solutions_and_sizes_by_seq(".json_pmE_MFE")[0]
    tab_all_F_MFE = get_nb_solutions_and_sizes_by_seq(".json_pmF_MFE")[0]
    data_to_plot = [tab_all_E_MEA, tab_all_F_MEA, tab_all_E_MFE, tab_all_F_MFE]

    all_data = []
    for i in range(len(data_to_plot)):
        for j in range(len(data_to_plot[i])):
            all_data.append(data_to_plot[i][j])

    min = np.amin(all_data)
    max = np.amax(all_data)
    median = np.median(all_data)

    fig = plt.figure()
    fig.set_size_inches(6, 3)
    ax = fig.add_axes([0, 0, 1, 1])

    labels = ['MEA + E', 'MEA + F', 'MFE + E', 'MFE + F']
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Répartition du nombre solutions pour chaque séquence du benchmark" + "\n minimum: " + str(min) + "    mediane: " + str(
            median) + "    maximum: " + str(max) + " (Pour l'ensemble des solutions)")
    ax.set_ylabel("nombre de solutions (structures secondaires prédites)")

    violins = ax.violinplot(data_to_plot, showmedians=True)

    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[partname]
        vp.set_edgecolor('blue')
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor("blue")
    plt.savefig('repartition_nb_solutions.png', dpi=200, bbox_inches='tight')

# Create a scatter plot showing the number of solutions according to the length of the sequence in the benchmark
def stats_nb_solutions_by_seq_length():
    list_name = get_list_name_bm()
    x = []
    y = []
    tab_all_E_MEA = get_nb_solutions_and_sizes_by_seq(".json_pmE_MEA")
    for i in range(len(tab_all_E_MEA[0])):
        x.append(tab_all_E_MEA[1][i])
        y.append(tab_all_E_MEA[0][i])
    tab_all_F_MEA = get_nb_solutions_and_sizes_by_seq(".json_pmF_MEA")
    for i in range(len(tab_all_F_MEA[0])):
        x.append(tab_all_F_MEA[1][i])
        y.append(tab_all_F_MEA[0][i])
    tab_all_E_MFE = get_nb_solutions_and_sizes_by_seq(".json_pmE_MFE")
    for i in range(len(tab_all_E_MFE[0])):
        x.append(tab_all_E_MFE[1][i])
        y.append(tab_all_E_MFE[0][i])
    tab_all_F_MFE = get_nb_solutions_and_sizes_by_seq(".json_pmF_MFE")
    for i in range(len(tab_all_F_MFE[0])):
        x.append(tab_all_F_MFE[1][i])
        y.append(tab_all_F_MFE[0][i])

    plt.scatter(x, y, s=50, c='blue', marker='o', edgecolors='black')
    plt.ylabel("nombre de solutions")
    plt.xlabel("longueur de la séquence")
    plt.title('nombre de structures prédites en fonction de la longueur de la séquence')
    plt.savefig('nb_solutions_en_fonction_seq.png', dpi=200, bbox_inches='tight')

"""stats_nb_motifs_in_result(".json_pmE_MEA")
stats_nb_motifs_in_result(".json_pmF_MEA")
stats_nb_motifs_in_result(".json_pmE_MFE")
stats_nb_motifs_in_result(".json_pmF_MFE")
stats_nb_motifs_all()
stats_nb_solutions()"""
stats_nb_solutions_by_seq_length()
