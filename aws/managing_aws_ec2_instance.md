# Managing your AWS EC2 instances


Here are a few helpful commands for starting and stopping your AWS instances when you come in/leave

First make sure AWS CLI is downloaded and on your ubuntu. You can follow the directions for LINUX found here [Installing or updating the latest version of the AWS CLI - AWS Command Line Interface](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) .

Once installed, update your credentials appropriately using aws configure.

I was having a few issues with my credentials so I just updated both my `~/.aws/credentials` and `~/.aws/config` files as such:

 

```
[default]
aws_access_key_id = ABCGDTYAHA
aws_secret_access_key = dy+-NAKUSDHFIOUASD
region = us-west-2
```

1.  Place the following code in a file named aws_start and then add the file to your `~/bin/` directory on your laptop. Make it executable: `chmod u+rwx aws_start`. Be sure to fill in the instance id with whatever one you use.

```
#!/bin/bash
aws ec2 start-instances  --instance-ids i-09aeiouperseph
```

2. Place the following code in a file named aws_stop and then add the file to your `~/bin/` directory on your laptop. Make it executable: `chmod u+rwx aws_stop` . Be sure to fill in the instance id with the same one as used in `aws_start`

```
#!/bin/bash
aws ec2 stop-instances  --instance-ids i-09aeioupersephone
```

Once done, open a new ubuntu terminal window and type aws_start. If it works you should see something like the following:

```
(base) pedrotorres@DESKTOP-BHKI6NK:~$ aws_start
{
    "StartingInstances": [
        {
            "CurrentState": {
                "Code": 0,
                "Name": "pending"
            },
            "InstanceId": "i-09aeioupersephone",
            "PreviousState": {
                "Code": 80,
                "Name": "stopped"
            }
        }
    ]
}
```

If you get the following error like I did:

`An error occurred (AuthFailure) when calling the StartInstances operation: AWS was not able to validate the provided access credentials`

the date on ubuntu might be off like mine was. Follow the instructions here on how to update the date on your ubuntu machine.

[Fix AWS CLI AuthFailure - Validate with access credentials](https://www.learnitguide.net/2018/04/fix-aws-cli-authfailure-validate-with.html) 

Once you have everything working you will be able to simply execute aws_start from your laptop when you come into work in the morning (or log onto your computer from home), and when you leave/are done you can execute aws_stop to terminate the instance while you are not using it. 

If you have multiple instances/machines that are stopped and want to use you can use the same process as above and update the names to reflect your machines such as aws_start_centrifuge and aws_stop_centrifuge ect..

## So here is a typical day:

1. execute `aws_start` from your ubuntu on windows environment 

2. give it ~30 seconds to boot up

3. ssh into the instance that you just started

4. Optional: if you need to mount_data you can do this within your instance

5. Do work

6. at the end of the day, execute `aws_stop` from your ubuntu on windows environment
