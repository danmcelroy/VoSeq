[![Build Status](https://travis-ci.org/carlosp420/VoSeq.svg)](https://travis-ci.org/carlosp420/VoSeq)

# Migrate VoSeq database
You need to dump your MySQL database into a XML file:

```shell
mysqldump --xml voseq_database > dump.xml
```

Then use our script to migrate all your VoSeq data into a PostGreSQL database.

```shell
python migrate_db.py dump.xml
```

It might issue a warning message:

```
WARNING:: Could not parse dateCreation properly.
WARNING:: Using empty as date for `time_edited` for code Your_Vocher_Code
```

It means that the creation time for your voucher was probably empty or similar
to `0000-00-00`. In that case the date of creation for your voucher will be
empty. This will not cause any trouble when running VoSeq. You can safely
ignore this message.
