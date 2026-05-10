// get the ninja-keys element
const ninja = document.querySelector('ninja-keys');

// add the home and posts menu items
ninja.data = [{
    id: "nav-about",
    title: "About",
    section: "Navigation",
    handler: () => {
      window.location.href = "/";
    },
  },{id: "nav-blog",
          title: "Blog",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/Blog/index.html";
          },
        },{id: "nav-publications",
          title: "Publications",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/Publications/";
          },
        },{id: "nav-cv",
          title: "CV",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/CV/";
          },
        },{id: "nav-photos",
          title: "Photos",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/photos/";
          },
        },{id: "nav-miscellaneous",
          title: "Miscellaneous",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/miscellaneous/";
          },
        },{id: "post-nematic-order-parameter-written-by-chatgpt-5-4-thinking",
        
          title: "Nematic Order Parameter (Written by ChatGPT 5.4 Thinking)",
        
        description: "A practical introduction to the nematic order parameter, including the Q tensor, the P2 Legendre polynomial, eigenvalue interpretation, and how to compute orientational order in practice.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2026/nematic/";
          
        },
      },{id: "post-openclaw-指北-韬小蛋的自述",
        
          title: "OpenClaw 指北：韬小蛋的自述",
        
        description: "",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2026/openclaw/";
          
        },
      },{id: "post-langevin-sampling-and-score-matching",
        
          title: "Langevin Sampling and Score Matching",
        
        description: "We often treat the &#39;normalization constant&#39; in deep learning as a mathematical headache. But by looking at it through the lens of physics, we find a solution. This post breaks down how Langevin dynamics allows us to bypass the intractable integral and generate data using the &#39;score function&#39;—physics&#39; answer to generative modeling.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/langevin-dynamic3/";
          
        },
      },{id: "post-practical-implementation-of-langevin-dynamics",
        
          title: "Practical Implementation of Langevin Dynamics",
        
        description: "An intuitive introduction to stochastic differential equations, starting from the Wiener process and progressing to Itô&#39;s lemma, the Fokker-Planck equation, and Langevin dynamics.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/langevin-dynamic2/";
          
        },
      },{id: "post-langevin-dynamics-from-wiener-process-to-stochastic-differential-equations",
        
          title: "Langevin Dynamics: From Wiener Process to Stochastic Differential Equations",
        
        description: "An intuitive introduction to stochastic differential equations, starting from the Wiener process and progressing to Itô&#39;s lemma, the Fokker-Planck equation, and Langevin dynamics.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/langevin-dynamics1/";
          
        },
      },{id: "post-opes-and-implementation-in-toy-model",
        
          title: "OPES and Implementation in Toy Model",
        
        description: "A practical guide to On-the-fly Probability Enhanced Sampling (OPES) with implementation details and a toy model demonstration.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/enhanced-sampling3/";
          
        },
      },{id: "post-collective-variables-and-free-energy-profile",
        
          title: "Collective Variables and Free Energy Profile",
        
        description: "Understanding collective variables (CVs) and how to estimate free energy profiles from molecular simulations using kernel density estimation.",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/enhanced-sampling2/";
          
        },
      },{id: "post-reweighting-technique-in-molecular-simulations-from-umbrella-sampling-to-mbar",
        
          title: "Reweighting Technique in Molecular Simulations: From Umbrella Sampling to MBAR",
        
        description: "An introduction to importance sampling and reweighting techniques in molecular simulations, covering umbrella sampling and MBAR (Multistate Bennett Acceptance Ratio).",
        section: "Posts",
        handler: () => {
          
            window.location.href = "/blog/2025/enhanced-sampling1/";
          
        },
      },{id: "books-dark-souls-remastered",
          title: 'DARK SOULS™: REMASTERED',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/darksoul/";
            },},{id: "books-hollow-knight-silksong-空洞骑士-丝之歌",
          title: 'Hollow Knight: Silksong / 空洞骑士：丝之歌',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/silksong/";
            },},{id: "books-dynasty-warriors-origins-真三国无双-起源",
          title: 'Dynasty Warriors: Origins / 真三国无双：起源',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/dynasty/";
            },},{id: "books-resident-evil-2-生化危机2重制版",
          title: 'Resident Evil 2 / 生化危机2重制版',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/re2/";
            },},{id: "books-resident-evil-requiem-生化危机-安魂曲",
          title: 'Resident Evil Requiem / 生化危机 安魂曲',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/re9/";
            },},{id: "books-resident-evil-7-biohazard-生化危机7",
          title: 'Resident Evil 7: Biohazard / 生化危机7',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/re7/";
            },},{id: "news-i-opened-my-personal-website-sparkles-smile",
          title: 'I opened my personal website! :sparkles: :smile:',
          description: "",
          section: "News",},{id: "projects-project-1",
          title: 'project 1',
          description: "with background image",
          section: "Projects",handler: () => {
              window.location.href = "/projects/1_project/";
            },},{id: "projects-project-2",
          title: 'project 2',
          description: "a project with a background image and giscus comments",
          section: "Projects",handler: () => {
              window.location.href = "/projects/2_project/";
            },},{id: "projects-project-3-with-very-long-name",
          title: 'project 3 with very long name',
          description: "a project that redirects to another website",
          section: "Projects",handler: () => {
              window.location.href = "/projects/3_project/";
            },},{id: "projects-project-4",
          title: 'project 4',
          description: "another without an image",
          section: "Projects",handler: () => {
              window.location.href = "/projects/4_project/";
            },},{id: "projects-project-5",
          title: 'project 5',
          description: "a project with a background image",
          section: "Projects",handler: () => {
              window.location.href = "/projects/5_project/";
            },},{id: "projects-project-6",
          title: 'project 6',
          description: "a project with no image",
          section: "Projects",handler: () => {
              window.location.href = "/projects/6_project/";
            },},{id: "projects-project-7",
          title: 'project 7',
          description: "with background image",
          section: "Projects",handler: () => {
              window.location.href = "/projects/7_project/";
            },},{id: "projects-project-8",
          title: 'project 8',
          description: "an other project with a background image and giscus comments",
          section: "Projects",handler: () => {
              window.location.href = "/projects/8_project/";
            },},{id: "projects-project-9",
          title: 'project 9',
          description: "another project with an image 🎉",
          section: "Projects",handler: () => {
              window.location.href = "/projects/9_project/";
            },},{
        id: 'social-email',
        title: 'email',
        section: 'Socials',
        handler: () => {
          window.open("mailto:%63%79%74%77%79%61%74%74@%6F%75%74%6C%6F%6F%6B.%63%6F%6D", "_blank");
        },
      },{
        id: 'social-github',
        title: 'GitHub',
        section: 'Socials',
        handler: () => {
          window.open("https://github.com/cytwyatt", "_blank");
        },
      },{
        id: 'social-scholar',
        title: 'Google Scholar',
        section: 'Socials',
        handler: () => {
          window.open("https://scholar.google.com/citations?user=WDXfwU0AAAAJ", "_blank");
        },
      },{
        id: 'social-steam_id',
        title: 'Steam_id',
        section: 'Socials',
        handler: () => {
          window.open("", "_blank");
        },
      },{
        id: 'social-wechat_qr',
        title: 'Wechat_qr',
        section: 'Socials',
        handler: () => {
          window.open("", "_blank");
        },
      },{
        id: 'social-zhihu_id',
        title: 'Zhihu_id',
        section: 'Socials',
        handler: () => {
          window.open("", "_blank");
        },
      },{
      id: 'light-theme',
      title: 'Change theme to light',
      description: 'Change the theme of the site to Light',
      section: 'Theme',
      handler: () => {
        setThemeSetting("light");
      },
    },
    {
      id: 'dark-theme',
      title: 'Change theme to dark',
      description: 'Change the theme of the site to Dark',
      section: 'Theme',
      handler: () => {
        setThemeSetting("dark");
      },
    },
    {
      id: 'system-theme',
      title: 'Use system default theme',
      description: 'Change the theme of the site to System Default',
      section: 'Theme',
      handler: () => {
        setThemeSetting("system");
      },
    },];
