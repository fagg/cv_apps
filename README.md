Work from Prof. Simon Lucey's 16-423 Designing Computer Vision Apps course.

https://github.com/slucey-cs-cmu-edu

Note:
- Drag opencv2.framework from old application to new
- Need to add Accelerate.framework from "General"->"Linked Frameworks and Libraries". Else it will complain about x86\_64 version.

To add Prof. Lucey's modules (http://stackoverflow.com/questions/4161022/git-how-to-track-untracked-content):
`git rm --cached AR_TennisBall/`
`rm -rf AR_TennisBall.git`
`git add AR_TennisBall/`
