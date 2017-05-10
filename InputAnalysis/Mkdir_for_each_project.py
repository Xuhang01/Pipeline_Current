import os
def main():
    Project_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Data/download/Project/"
    Output_dir = "/f/hang/Data/GeneFusion/Database/ICGCdataset/Project_output/"
    for Project_ID in os.listdir(Project_dir):
        if not os.path.exists(os.path.join(Output_dir, Project_ID)):
            os.mkdir(os.path.join(Output_dir, Project_ID))
if __name__ == '__main__':
    main()