try:
    from my_credentials import username, pw_local, pw_remote
except:
    print(
        '''Create a file "my_credentials.py" containing:
        username = 'yourname'
        pw_local = 'your local password'
        pw_remove = 'your remote password'
        ''')
    raise


from core_tools.data.SQL.connect import (
        set_up_local_storage,
        set_up_remote_storage,
        set_up_local_and_remote_storage
        )

def setup():
    set_up_local_storage(
            username, pw_local, "sds_test",
            "test_Tektronix", "test_system", "no"
            )

def setup_local_and_remote():
    set_up_local_and_remote_storage(
            'vanvliet.qutech.tudelft.nl', 5432,
            username, pw_local, "sds_test",
            username, pw_remote, "sds_test",
            "test_Tektronix", "test_system", "no"
            )

def setup_local_other(db_name, project=None, setup=None, sample=None):
    set_up_local_storage(
            username, pw_local, db_name,
            project, setup, sample
            )

def setup_local_and_remote_other(db_name, project=None, setup=None, sample=None,
                                 local_readonly=False, remote_readonly=True):
    set_up_local_and_remote_storage(
            'vanvliet.qutech.tudelft.nl', 5432,
            username, pw_local, db_name,
            username, pw_remote, db_name,
            project, setup, sample,
            local_readonly, remote_readonly
            )

def setup_remote_other(db_name, project=None, setup=None, sample=None, readonly=True):
    set_up_remote_storage(
            'vanvliet.qutech.tudelft.nl', 5432,
            username, pw_remote, db_name,
            project, setup, sample,
            readonly
            )
