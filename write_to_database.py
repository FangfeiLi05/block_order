import psycopg2
import shutil
import os
import argparse


def write_to_database_step1(filename):
    """
    """
    print('Writing to database...')
    filename_tmp = '/tmp/{}'.format(filename.split('/')[-1])
    shutil.move(filename, filename_tmp)
    change_right = 'chmod 777 {}'.format(filename_tmp)
    os.system(change_right)

    conn = psycopg2.connect(database='my_database',
                            user='postgres',
                            password='870525',
                            host='localhost',
                            port='5432')
                            
    #conn = psycopg2.connect(database='our_database',
    #                        user='postgres',
    #                        password='pro3165',
    #                        host='localhost',
    #                        port='5432')

    cur = conn.cursor()
    cur.execute('CALL write_to_blockOrder_block_adapter (%s);', (filename_tmp,))
    cur.execute('COMMIT')
    cur.close()
    conn.close()
    
    shutil.move(filename_tmp, filename)
    print('Done!')
        


def main():
    """
    """
    parser = argparse.ArgumentParser(description='Script for writing to databse',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-s1', '--step1', type=str, required=True,
                        help='a .csv file')
                    
    args = parser.parse_args()
    
    print('Writing to database...')
    filename_step1 = args.step1
    write_to_database_step1(filename_step1)
    print('Done!')
    
if __name__ == "__main__":
    main()
