import time
from matplotlib import pyplot as plt


BLOSUM62 = {}
amino = "ARNDCQEGHILKMFPSTWYV"    # All 20 Amino Acids representing the keys for the BLOSUM62 Dictionary 

matrix =[[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],

 [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],

 [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],

 [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],

 [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],

 [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],

 [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],

 [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],

 [-2, 0, 1,-1,-3,-0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],

 [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],

 [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],

 [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],

 [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],

 [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],

 [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],

 [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],

 [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],

 [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],

 [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],

 [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]
]

# Function for the initialization of BLOSUM62 dictionary
def create_blosum():
    for i in range(len(amino)):
        for j in range(len(amino)):
            temp = (amino[i],amino[j])
            BLOSUM62[temp] = matrix[i][j]
    return

# Needleman-Wunsch algorithm with affine gap penalties for global alignment of sequences
def needleman_wunsch_affine(s1, s2, gap_open, gap_extend):
    n = len(s1)
    m = len(s2)
    d = [[0 for i in range(m+1)] for j in range(n+1)]       # Table d indicates the matches/mismatches
    a = [[0 for i in range(m+1)] for j in range(n+1)]       # Table a indicates the gaps extended in s1
    b = [[0 for i in range(m+1)] for j in range(n+1)]       # Table b indicates the gaps extended in s2
    t = [[0 for i in range(m+1)] for j in range(n+1)]       # Table t is utilized for tracing back the resulting alignment
    string_show = False    # If True, the resultant strings will be printed in the terminal

    # Initializing d matrix
    d[0][0] = 0
    d[0][1] = gap_open + gap_extend
    d[1][0] = gap_open + gap_extend
    for i in range(2,len(d)):
        d[i][0] = d[i-1][0] + gap_extend
    for i in range(2,len(d[0])):
        d[0][i] = d[0][i-1] + gap_extend

    # Inititalizing a matrix
    a[0][1] = gap_open + gap_extend
    for i in range(len(a)):
        a[i][0] = -10000000000
    for i in range(2,len(a[0])):
        a[0][i] = a[0][i-1] + gap_extend
    
    # Initializing b matrix
    b[1][0] = gap_open + gap_extend
    for i in range(len(b[0])):
        b[0][i] = -10000000000
    for i in range(2,len(b)):
        b[i][0] = b[i-1][0] + gap_extend

    # Initializing t matrix
    for i in range(len(t)):
        t[i][0] = 'b'
    for i in range(len(t[0])):
        t[0][i] = 'a'
    t[0][0] = '-'

    # Updating all the dp (dynamic programming) matrices according to the BLOSUM62 matrix
    for i in range(1,len(d)):
        for j in range(1,len(d[0])):
            b[i][j] = max(b[i-1][j] + gap_extend, d[i-1][j] + gap_extend + gap_open)
            a[i][j] = max(a[i][j-1] + gap_extend, d[i][j-1] + gap_extend + gap_open)
            scr = BLOSUM62[(s1[i-1],s2[j-1])]
            d[i][j] = max(d[i-1][j-1] + scr, a[i][j], b[i][j])
            if d[i][j] == d[i-1][j-1] + scr:
                t[i][j] = 'd'
            elif d[i][j] == a[i][j]:
                t[i][j] = 'a'
            elif d[i][j] == b[i][j]:
                t[i][j] = 'b'

    # Initializing the resultant strings 
    res1 = ""
    res2 = ""
    i = n
    j = m

    # Trace back to obtain the resultant strings
    while(i > 0 or j > 0):
        if t[i][j] == 'd':
            res1 += s1[i-1]
            res2 += s2[j-1]
            i -= 1
            j -= 1
        elif t[i][j] == 'b':
            res1 += s1[i-1]
            res2 += '_'
            i -= 1
        elif t[i][j] == 'a':
            res1 += '_'
            res2 += s2[j-1]
            j -= 1

    # Reversing the resultant strings
    res1 = res1[::-1]
    res2 = res2[::-1]

    # Outputting the resulting strings if string_show is True
    if(string_show):
        print("\nResultant strings: \n")
        print(res1)
        print(res2)

    # Return the resultant strings res1 and res2, and the score of the alignment given by d[n][m]
    return res1, res2, d[n][m]

# Function to calculate the identity of the calculated optimal alignment
def show_identity(s1,s2):
    iden = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            iden += 1
    #print("Identity: " + str(round(float(iden/len(s1)),3)))
    return float(iden/len(s1))

# Function to calculate the similarity of the calculated optimal alignment
def show_similarity(s1,s2):
    simi = 0
    for i in range(len(s1)):
        if s1[i] == "_" or s2[i] == "_":
            continue
        if BLOSUM62[(s1[i],s2[i])] > 0:
            simi += 1
    #print("Similarity: " + str(round(float(simi/len(s1)),3)))
    return float(simi/len(s1))

if __name__ == "__main__":

    show_time = True

    start_time = time.time()

    # Call for function to initialize the BLOSUM62 matrix
    create_blosum()

    # Input the strings
    s1 = input()    
    s2 = input()

    # Initializing variables and the 2-D arrays for storing the respective values for all the gap penalty iterations
    iden = [[0 for x in range(2000)] for y in range(2000)]
    sim = [[0 for x in range(2000)] for y in range(2000)]
    scores = [[0 for x in range(2000)] for y in range(2000)]

    max_score = -1000000000
    min_score = 1000000000

    max_sim = -1000000000
    min_sim = 1000000000

    max_iden = -1000000000
    min_iden = 1000000000

    counter1 = 0    # Used for indexing
    counter2 = 0

    # The lower bound for (gap_open / gap_extend) in order to ensure that the gap_open is significantly greater than gap_extend
    kappa = 5.0

    # Calculating the values for every gap_open and gap_extend penalties
    for index1 in range(-200,0,1): 
        counter2 = 0    

        i = (float)(index1/10)    # Increments are done by 1/10 = 0.1

        for index2 in range(index1,0,1):    # Since (kappa.|gap_extend|) < |gap_open|, the for loop is modified accordingly

            j = (float)(index2/(10*kappa))    # Increments are done by 1/(10*kappa)

            res1, res2, score = needleman_wunsch_affine(s1,s2,i,j)    # Running the Needleman-Wunsch algorithm for every iteration

            # Updating the matrices with the values obtained
            scores[counter1][counter2] = score
            iden[counter1][counter2] = show_identity(res1,res2)
            sim[counter1][counter2] = show_similarity(res1,res2)

            # Finding the maximum and minimum values obtained for each variable which will be used later for normalization
            max_score = max(max_score,score)
            min_score = min(min_score,score)

            max_sim = max(max_sim,sim[counter1][counter2])
            min_sim = min(min_sim,sim[counter1][counter2])

            max_iden = max(max_iden,iden[counter1][counter2])
            min_iden = min(min_iden,iden[counter1][counter2])

            counter2 += 1
        counter1 += 1

    print("Gap open and Gap extend caluclated for every iteration")

    # Iteration to normalize all the variables to take values from 0 to 1 by utilising the maximum and the minimum values of each
    for counter1 in range(201):

        for counter2 in range(201):

            if max_score == min_score:    # If the values obtained are the same throughout, we choose a different normalization constant
                scores[counter1][counter2] = 0
            else:    # Else, normalize the values to constrain them between 0 and 1 (inclusive)
                scores[counter1][counter2] = (float)((scores[counter1][counter2] - min_score)/(max_score-min_score))

            if max_sim == min_sim:
                sim[counter1][counter2] = max_sim
            else:
                sim[counter1][counter2] = (float)((sim[counter1][counter2] - min_sim)/(max_sim - min_sim))

            if max_iden == min_iden:
                iden[counter1][counter2] = max_iden
            else:
                iden[counter1][counter2] = (float)((iden[counter1][counter2] - min_iden)/(max_iden - min_iden))

    # Initializing the 2-D array which will be used to store new scores which are calculated using the new scoring system
    new_scores = [[0 for x in range(202)] for y in range(202)]

    # Initializing the indices
    counter1 = 0
    counter2 = 0

    # Calculate the new scores by using the defined new scoring system and iterating over all the values of the parameters
    # Iterating over all the values of the parameters with an increment of 0.01
    for gamma_temp in range(101):   

        counter2 = 0
        gamma = (float)(gamma_temp/100)

        for beta_temp in range(101-gamma_temp):

            beta = (float)(beta_temp/100)
            alpha = 1 - gamma - beta    # Definining the three parameters by iterating over all values

            if alpha < 0 or beta < 0 or gamma < 0 or alpha > 1 or beta > 1 or gamma > 1:
                continue

            for i in range(201):

                for j in range(201): 

                    # Update and add all the new scores for every pair of gap penalties to update the new scoring matrix according to the scoring system
                    new_scores[counter1][counter2] += (float)(alpha * iden[i][j] + beta * sim[i][j] + gamma * scores[i][j])

            # Calculate the mean of the new scores for every value of alpha, beta, and gamma (and store it in the new_scores matrix)                   
            new_scores[counter1][counter2]/=(40000)

            counter2 += 1

        counter1 += 1

    print("Calculated new scores")

    # Initialize the maximal/ best new score to be calculated
    best_new_score = -10000000

    # Finding the maximal new score and thus, computing the corresponding alpha, beta, and gamma values for that maximal new score
    for i in range(201):

        for j in range(201):

            if best_new_score < new_scores[i][j]: 

                best_new_score = new_scores[i][j]

                # Calculating the final values of parameters to be used
                final_gamma = (float)(i/200)
                final_beta = (float)(j/200)
                final_alpha = (1 - final_gamma - final_beta)

    show_all_parameters = False    # If all the viable values for parameters are to be printed, set as True. Otherwise, set as False

    print("Calculated the new parameters to be used")

    # Calculating the final values of parameters that can be used alternatively with the same effect on the scoring
    for i in range(201):

        for j in range(201):

            if best_new_score == new_scores[i][j]:

                final_temp_gamma = (float)(i/200)
                final_temp_beta = (float)(j/200)
                final_temp_alpha = (1 - final_temp_gamma - final_temp_beta)

                if(show_all_parameters):
                    print("All the viable values of parameters:")
                    print("Alpha: "+str(final_temp_alpha) + ", Beta: " + str(final_temp_beta) + ", Gamma: "+str(final_temp_gamma))

    # Initialization of the final gap penalties to be outputted, and the maximal score to be calculated using the final parameter values
    maximum_recalculated_score = -100000000
    best_gap_open = -20.0
    best_gap_extend = -20.0

    for index1 in range(-200,0,1): 

        counter2 = 0
        i = (float)(index1/10)

        for index2 in range(index1,0,1):

            j = (float)(index2/(10*kappa))

            # Calculate the new scores according to the final calculated values of the parameters
            val = (final_alpha * iden[counter1][counter2] + final_beta * sim[counter1][counter2] + final_gamma * scores[counter1][counter2])

            # Store the best gap penalties to be outputted for the corresponding maximum recalculated score
            if val > maximum_recalculated_score:

                maximum_recalculated_score = val
                best_gap_open = i
                best_gap_extend = j

            counter2 += 1

        counter1+=1

    #Printing the final output
    print("Final Output:")
    print(final_alpha,final_beta,final_gamma)
    print(best_gap_open,best_gap_extend)

    final_s1, final_s2, final_score = needleman_wunsch_affine(s1,s2,best_gap_open,best_gap_extend)
    final_simi = show_similarity(final_s1,final_s2)
    final_iden = show_identity(final_s1,final_s2)

    new_final_score = final_alpha*final_iden + final_beta*final_simi + final_gamma*final_score

    print(final_s1)
    print(final_s2)
    print(final_iden, final_simi, final_score, new_final_score)

    end_time = time.time()

    # Print the time elapsed, if show_time == True
    if(show_time):
        print("Time Elapsed :" + str(end_time - start_time) + " sec")