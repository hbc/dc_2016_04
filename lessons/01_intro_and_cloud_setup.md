---
layout: page
title: "The Shell"
comments: true
date: 2014-07-30
---

# The Shell

## What is the shell?

The *shell* is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination. **The *shell* is an interpreter that helps translate our input into computer language.**

There are many reasons to learn about the shell.

* For most bioinformatics tools, you have to use the shell. There is no
graphical interface. If you want to work in metagenomics or genomics you're
going to need to use the shell.
* The shell gives you *power*. The command line gives you the power to do your work more efficiently and
more quickly.  When you need to do things tens to hundreds of times,
knowing how to use the shell is transformative.
* Computational resources that can handle large datasets, such as remote computers or cloud computing, require a working knowledge of Unix.


![Automation](../img/gvng.jpg)

  Unix is user-friendly. It's just very selective about who its friends are.

Today we're going to go through how to access Unix/Linux and some of the basic
shell commands.


##Setting up

For the duration of this workshop, we will be accessing the shell on a remote computer using the **Amazon's Elastic Compute Cloud (EC2)**. Remote machines are invaluable to scientists working with large datasets when...

- Their computer does not have enough resources to run the desired analysis (memory, processors, disk space, network bandwidth).
- Their computer is taking hours or days to get through an analysis.
- They cannot install software on their computer (application does not have support for their operating system, conflicts with other existing applications, etc.)

We will spend most of our time learning about the basics of the shell by manipulating some experimental data. Since we are going to be working with this data on the cloud, we first need to log onto a remote computer.

### Connecting to a remote computer using Amazon EC2


This is the first and last place in these lessons where it will matter if you are using PC, Mac, or Linux. After we connect, we will all be on the same operating system/computing environment. 

> To save time, your instructor will have launched an remote computer (instance) for you prior to the workshop. If you are following these lessons on your own, or after the workshop see the lesson on [cloud computing](https://github.com/datacarpentry/cloud-genomics/tree/gh-pages/lessons) for instructions on how to do this yourself. 

#####**User Credentials**
Credentials are case sensitive:

- Username: dcuser
- Password: data4Carp

#####**Connecting using PC**
*Prerequisites*: You must have an SSH client. There are several free options and we will use PuTTY [[Download Putty.exe](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)]

1. Open PuTTY; In the 'Host Name (or IP address)' section paste in the IP address provided by your instructor (or the ip address of an instance you have provisioned yourself). *Keep the default selection 'SSH' and Port (22)*. <br>
![Putty Image](../img/putty_screenshot_1.png)
2. Click 'Open' and you will be presented with a security warning. Select 'Yes' to continue to connect. <br>
![Putty security screen](../img/putty_screenshot_2.png)
3. In the final step, you will be asked to provide a login and password. **Note:** When typing your password, it is common in Unix/Linux not see see any asterisks (e.g. ****) or moving cursors. Just continue typing.<br> 
![Putty login](../img/putty_screenshot_3.png)
4. You should now be connected!

#####**Connecting using Mac/Linux**
*Prerequisites*: Mac and Linux operating systems will already have terminals installed. Simply search for 'Terminal' and/or look for the terminal icon.<br> 
![terminal icon](../img/terminal.png)


1. open the terminal and type the following command substituting 'ip_address' for the ip address your instructor will provide (or the ip address of an instance you have provisioned yourself). *Be sure to pay attention to capitalization and spaces*

        $ ssh dcuser@ip_address
        
2. You will receive a security message that looks something like the message below. Type 'yes' to proceed.

        The authenticity of host 'ec2-52-91-14-206.compute-1.amazonaws.com (52.91.14.206)' can't be established. ECDSA key fingerprint is SHA256:S2mMV8mCThjJHm0sUmK2iOE5DBqs8HiJr6pL3x/XxkI. Are you sure you want to continue connecting (yes/no)?

3. In the final step, you will be asked to provide a login and password. **Note:** When typing your password, it is common in Unix/Linux not see see any asterisks (e.g. ****) or moving cursors. Just continue typing.
4. You should now be connected!



###**Verifying your connection and environment** 

When you connect, it is typical to receive a welcome screen. The Data Carpentry Amazon instances display this message upon connecting:


```
Welcome to Ubuntu 14.04.3 LTS (GNU/Linux 3.13.0-48-generic x86_64)

 * Documentation:  https://help.ubuntu.com/

  System information as of Sun Jan 24 21:38:35 UTC 2016

  System load:  0.0                Processes:           151
  Usage of /:   48.4% of 98.30GB   Users logged in:     0
  Memory usage: 6%                 IP address for eth0: 172.31.62.209
  Swap usage:   0%

  Graph this data and manage this system at:
    https://landscape.canonical.com/

  Get cloud support with Ubuntu Advantage Cloud Guest:
    http://www.ubuntu.com/business/services/cloud

12 packages can be updated.
10 updates are security updates.


Last login: Sun Jan 24 21:38:36 2016 from
```

You should also have a blinking cursor awaiting your command


```bash
dcuser@ip-172-31-62-209 ~ $
```
