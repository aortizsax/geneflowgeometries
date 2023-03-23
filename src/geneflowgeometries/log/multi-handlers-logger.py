import logging
from logging.handlers import RotatingFileHandler


def setup_logger(config_dict):
    MAX_BYTES = 10000000  # Maximum size for a log file
    BACKUP_COUNT = 9  # Maximum number of old log files

    # The name should be unique, so you can get in in other places
    # by calling `logger = logging.getLogger('com.dvnguyen.logger.example')
    logger = logging.getLogger("geneflowgeometries.log")
    logger.setLevel(
        logging.INFO
    )  # the level should be the lowest level set in handlers

    log_format = logging.Formatter("[%(levelname)s] %(asctime)s - %(message)s")

    stream_handler = logging.StreamHandler(host="localhost")
    stream_handler.setFormatter(log_format)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    info_handler = RotatingFileHandler(
        "info.log", maxBytes=MAX_BYTES, backupCount=BACKUP_COUNT
    )
    info_handler.setFormatter(log_format)
    info_handler.setLevel(logging.INFO)
    logger.addHandler(info_handler)

    error_handler = RotatingFileHandler(
        "error.log", maxBytes=MAX_BYTES, backupCount=BACKUP_COUNT
    )
    error_handler.setFormatter(log_format)
    error_handler.setLevel(logging.ERROR)
    logger.addHandler(error_handler)


# if __name__ == '__main__':
#  setup_logger()
#  for i in range(0, 1000):
#    logger.info('This is a message {}'.format(i))
#    if i % 5 == 0:
#      logger.error('THis is a error {}'.format(i))


# use socket handler for saving to local host
# change to one column at the end of the experiment
#
# %(asctime)s %(levelname)s: %(message)s
# logging.info('LOGISTICSseed,PARAMETERS...., ANALYSISDATA...')per snapshot per trial
#


# experiment log
# logging.basicConfig(filename='experiment.log',
#                   filemode='w',
#                  level=logging.DEBUG,
#                 format='%(asctime)s %(levelname)s: %(message)s')


# use for confinguring log
# import logging
# import logging.config
# logging.config.fileConfig('logging.conf')


# .local log
# logging.basicConfig(filename='main.log',
#                   encoding='utf-8',
#                  level=logging.DEBUG,
#                 format='%(asctime)s %(levelname)s: %(message)s')
# logging.warning('is when this event was logged.')
# logging.debug('This message should go to the log file')
# logging.info('So should this')
# logging.warning('And this, too')#use for migrate value high
# logging.error('And non-ASCII stuff, too, like Øresund and Malmö')


# use in simulators
## mylib.py
# import logging
#
# def do_something():
#    logging.info('Doing something')
#
##logging.warning('%s before you %s', 'Look', 'leap!')


# now = str(datetime.datetime.now())

# print(now.split(" ")[0])
# date = now.split(" ")[0]
# main_log_row = now.split(" ")[0] + ","
###############################
#############################
# main_log_row += simulate_what + ","
# main_log_row += Geometry + ","
# main_log_row += str(number_of_chromosomes) + ","
# main_log_row += str(number_of_ploidy) + ","
# main_log_row += str(number_of_demes) + ","
# main_log_row += str(migration_rate) + ","
# main_log_row += str(mutation_rate) + ","
# main_log_row += str(sequence_length) + ","
# main_log_row += str(simT) + ","
# main_log_row += str(snapshot_times) + ","
######################

# f = open("main.log", "a")
# f.write(main_log_row)
# f.close()
