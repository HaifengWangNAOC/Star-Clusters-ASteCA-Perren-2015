

This [link](https://github.com/Gabriel-p/asteca/compare/791a238...cfbe6b9)
contains all the commits done to generalize the photometric
systems used, made between 09/10 and 04/11.
All this commits are applied inside this branch.

The following commits have been applied into the `master` branch
either using `git cherry-pick` or manually (from oldest to newest):

* [35d6e86](https://github.com/Gabriel-p/asteca/commit/35d6e86)
* [5b1ea69](https://github.com/Gabriel-p/asteca/commit/5b1ea69)
* [2a74b45](https://github.com/Gabriel-p/asteca/commit/2a74b45)
* [5f5976c](https://github.com/Gabriel-p/asteca/commit/5f5976c)
* [242b211](https://github.com/Gabriel-p/asteca/commit/242b211)
* [8dc9e08](https://github.com/Gabriel-p/asteca/commit/8dc9e08)
* [26eabba](https://github.com/Gabriel-p/asteca/commit/26eabba)
* [548f5db](https://github.com/Gabriel-p/asteca/commit/548f5db)
* [43cf465](https://github.com/Gabriel-p/asteca/commit/43cf465)

Once v0.1.1 is stable and released, create a new branch
and add one by one all the changes made in the commits
shown in the above link, until the new branch is as simlar
as this one as possible.

The command:

    git difftool branch1..branch2

can be used to check the differences between branches.

This branch (the old one) can be deleted after the above step
is done, making sure all the changes were moved into the
new branch.

After that, keep working (now on the new branch) until
the generalization is done. in the meantime, keep
the code in the master branch separated.

When the generalization is done, merge into master and
publish the code as v0.2.0.
