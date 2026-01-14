# Windows 下 OpenSSH 密钥权限修复与多跳配置标准

> **适用场景**: Windows 10/11 作为客户端，通过跳板机（如 Tailscale/Bastion）连接内网目标服务器。解决著名的 `WARNING: UNPROTECTED PRIVATE KEY FILE!` 报错。

## 核心痛点

OpenSSH 对私钥文件的权限要求极其严格（等同于 Linux 的 `chmod 600`）。但 Windows 的 NTFS 文件系统权限默认继承自父目录，导致私钥文件往往拥有“Authenticated Users”等组的读取权限，从而被 SSH 客户端拒绝加载。

本方案提供一套标准的 PowerShell 操作流，一键修复权限并配置 `ProxyJump`。

---

## 第一式：私钥权限修复 (The "chmod 600" on Windows)

在 PowerShell 中执行。将 `<KeyName>` 替换为你的私钥文件名（如 `id_rsa_jump`）。

```powershell
# 定义私钥路径 (请修改此处)
$KeyPath = "$env:USERPROFILE\.ssh\<KeyName>"

# 1. 【重置】: 清除所有显式权限，回归继承状态
icacls $KeyPath /reset

# 2. 【授权】: 仅授予当前登录用户“读取(R)”权限 (相当于 chmod 400/600 的所有者部分)
icacls $KeyPath /grant:r "$($env:USERNAME):(R)"

# 3. 【隔离】: 移除继承权限 (关键步骤！彻底切断 System/Administrators 组的默认访问，模拟 Linux 严格权限)
icacls $KeyPath /inheritance:r
```

> **验证**: 执行完后，该文件在“属性 -> 安全”中应只显示当前用户一个条目。

---

## 第二式：密钥托管 (SSH-Agent)

为了避免每次连接都输入密码，以及为了支持 `ForwardAgent`，建议启用 Windows 自带的 SSH Agent 服务。

```powershell
# 1. 设置服务自启动 (防止重启失效)
Set-Service -Name ssh-agent -StartupType Automatic

# 2. 启动服务
Start-Service ssh-agent

# 3. 加载私钥 (只需执行一次)
ssh-add "$env:USERPROFILE\.ssh\<KeyName>"
```

---

## 第三式：多跳路由配置 (ProxyJump)

编辑用户目录下的 SSH 配置文件：`~/.ssh/config`。

### 场景拓扑
`Local (Windows)` -> `JumpHost (跳板机)` -> `Target (内网目标)`

```ssh
# --- 跳板机配置 ---
Host jump-server
    HostName jump.example.com          # 跳板机公网 IP 或域名
    User root                          # 跳板机用户名
    Port 22                            # 跳板机端口
    IdentityFile "~/.ssh/id_rsa_jump"  # 跳板机私钥
    AddKeysToAgent yes                 # 自动加入 Agent
    # ForwardAgent yes                 # [可选] 仅当你需要在跳板机上再用 ssh 命令连其他机器时开启 (注意安全风险)
    ServerAliveInterval 60             # 防止空闲断连

# --- 目标机配置 (核心魔法) ---
Host target-server
    HostName 192.168.x.x               # 目标机内网 IP
    User target_user                   # 目标机用户名
    Port 22
    
    # [核心指令] 指定通过 jump-server 跳转
    # 本地 ssh 会先建立到 jump-server 的隧道，再通过隧道直连目标
    ProxyJump jump-server
    
    # 指定目标机的私钥 (注意：该私钥位于【本地 Windows】，而非跳板机上)
    IdentityFile "~/.ssh/id_rsa_target" 
```

### 部署目标机公钥的小技巧

如果你无法直接连接目标机，如何把公钥放上去？通过跳板机“隔山打牛”：

```powershell
# 语法: ssh <跳板机> "命令"
# 作用: 在目标机上创建目录并写入公钥 (需手动替换 <PublicKeyContent>)
ssh jump-server "ssh target_user@192.168.x.x 'mkdir -p ~/.ssh && chmod 700 ~/.ssh && echo <PublicKeyContent> >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys'"
```
