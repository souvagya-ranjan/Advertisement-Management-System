from collections import defaultdict
def compare_text_files(file1, file2):

    mapp1 = defaultdict(list)
    mapp2 = defaultdict(list)

    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        # Read the first line (integer n) from both files
        n1 = int(f1.readline().strip())
        n2 = int(f2.readline().strip())

        print(n1, n2)
        # Compare the values of n from both files
        if n1 != n2:
            print("n values from both files are not equal")
            return False

        # Read the next 2*n lines from both files and compare line by line
        for _ in range(n1):
            # print(_)
            v1 = int(f1.readline().strip())
            v2 = int(f2.readline().strip())

            # Compare the lines based on the criteria you provided
            # if int(v1) != int(v2):
            #     return False
            # print(v1,v2)
            ele1 = f1.readline().strip()
            ele2 = f2.readline().strip()
            
            arr1 = list(map(int,ele1.split()))
            arr2 = list(map(int,ele2.split()))
            arr1.sort()
            arr2.sort()


            # if arr1 != arr2:
            #     return False

            mapp1[v1] = arr1
            mapp2[v2] = arr2

            
        keys1 = sorted(mapp1.keys(), key = lambda x: len(mapp1[x]))
        keys2 = sorted(mapp2.keys(), key = lambda x: len(mapp2[x]))
        # if keys1!= keys2:
        #     print("keys set are not equal")
        #     return False
        counter = 0
        # print("keys1",keys1)
        # print("keys2",keys2)
        for i in range(n1):
            key = keys1[i]
            # counter += 1
            if mapp1[key] != mapp2[key]:
                print("counter = ", i)
                print(key)
                print(len(mapp1[key]), mapp1[key])
                print(len(mapp2[key]), mapp2[key])
                return False
        # If all lines are identical, return True
    return True

# Usage example:
file1 = "A3_test/test3/task2_output3_verbose.txt"  # Replace with the actual filename of file1
file2 = "Output/ayushT2_output3.txt"  # Replace with the actual filename of file2

if compare_text_files(file1, file2):
    print("Both files are identical.")
else:
    print("Files are not identical.")
