https://www.googleapis.com/auth/drive.readonly in https://developers.google.com/oauthplayground/?code=4/0Adeu5BVkQCoY8d7RJDPbpNQTVSvev-x2Tvm-yrV0pY3HF8tEiejBoSBgFh-zqFncl2RT7Q&scope=https://www.googleapis.com/auth/drive.readonly

cd /home/baranov_lab/10X/NF1/Aboozar_new

curl -H "Authorization: Bearer ya29.a0AfB_byC1mek9m4YyoF5bAKd0Qcv3t4I14lE_CCQTDoWs1dZYok9ComUtLr-rVBHe3_OVK2c9pRix_8s9EeC-sAbI4p4Kkz_8tG2mz9o0V_rhrc3QuboaPNCBqtHkjEptETkUAJ8uhMFRIO_5uJb9qEECyivzFbKdoHa1DAaCgYKATUSARMSFQHsvYlsN2KT1YY358n7bzvbggJudg0173" https://www.googleapis.com/drive/v3/files/1CJBVpN7I5FdJ4-x41-XXnH0P90HPIGhu/?alt=media -o bam9.bam
sudo samtools sort -t CB -O BAM -o/home/baranov_lab/10X/NF1/Aboozar_new/4/4/outs/cellsorted_possorted_genome_bam.bam /home/baranov_lab/10X/NF1/Aboozar_new/4/4/outs/possorted_genome_bam.bam
velocyto run10x -m /home/baranov_lab/10X/mus_musculus/mm10_rmsk.gtf /home/baranov_lab/10X/NF1/Aboozar_new/6/6 /home/baranov_lab/10X/mus_musculus/gencode.vM27.annotation.gtf

#make sure to reestablish connection every hour - EK
