# pp_LU_Decomposition

Use C and MPI to generate a random large square matrix and find the LU decomposition of a large square matrix.


**RUN:**

**_.pbs_** files are the job files that are submitted to the cluster where we specify the number of nodes, name of the executable file, and all other configurations for the job to run.

**_.pl_** files are the script files that automate the whole process by executing the c program and submitting the .pbs job to the cluster.

**_.err_** files give any error after executing the job.

**_.out_** files are the output files after the job being executed.

Run the .pl file to start the program execution on the server.
