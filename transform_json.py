import json
import os



if __name__=="__main__":

    ##json import
    dir = os.getcwd() + "/data/modules/CaRNAval/"
    try:
        file = open(dir+"dataset.json")

    except OSError:
        print("File not found : " + dir + "dataset.json")

    else:
        obj = json.load(file)
        file.close()
        print("Dataset loaded")


        ##Creation of a file for each RIN
        try:
            os.mkdir(dir + "Subfiles")
        except OSError:
            print ("Creation of the directory %s failed" % (dir + "Subfiles") + " : maybe it already exists ?")
        else:
            print ("Successfully created the directory %s " % (dir + "Subfiles"))

        header = "pos;k;seq\n"

        for motif in obj:
            f = open( dir + "Subfiles/" + str(motif["uid"]) + ".txt"  ,  "w+" )
            f.write(header)

            num_nt = 0

            for component in motif["l_graphs"]:
                k = len(component["nodes"])
                pos = (num_nt, num_nt+k-1)
                num_nt += k

                seq = ""

                for node in component["nodes"]:
                    seq += node["real_nt"]

                data_comp = str(pos[0]) + "," + str(pos[1]) + ";" + str(k) + ";" + str(seq) + "\n"

                f.write(data_comp)

            f.close()

        print("Successfully parsed "+dir+"dataset.json")

